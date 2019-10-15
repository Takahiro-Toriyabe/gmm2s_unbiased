/* Main interface */
capture program drop gmm2s_unbiased
program define gmm2s_unbiased, eclass
	syntax varlist(min=1 max=1 numeric) [aw] [if] [in], ///
		treatment(varname) /// Endogeneous variable
		iv_list(varlist min=1 numeric) /// Set of instruments
		[controls(varlist numeric fv ts)] /// Control variables
		[nsim(integer 1000000)] /// #simulation in calculation of BR estimate
		[vce(string)] /// VCE type
	
	marksample touse

	/* Partial out y, d, z */
	tempvar u_y u_d
	qui _partial_out `varlist' `controls' [`weight'`exp'] ///
		if `touse', gen(`u_y')
	qui _partial_out `treatment' `controls' [`weight'`exp'] ///
		if `touse', gen(`u_d')
	
	tempname cnt iv_list_u
	local `cnt' = 0
	foreach var in `iv_list' {
		tempvar u_z`++`cnt''
		qui _partial_out `var' `controls' [`weight'`exp'] ///
			if `touse', gen(`u_z``cnt''')
		local `iv_list_u' "``iv_list_u'' `u_z``cnt'''"
	}

	/* Apply unbiased GMM 2SLS estimation */
	tempname opt
	local `opt' ""
	if "`nsim'" != "" & ``cnt'' != 1 {
		local `opt' "``opt'' nsim(`nsim')"
	}
	if "`vce'" != "" {
		local `opt' "``opt'' vce(`vce')"
	}
	
	tempname command
	local `command' = "iv_unbiased" + (``cnt'' != 1) * "_gmm"
	qui ``command'' `u_y' [`weight'`exp'] if `touse', ///
		treatment(`u_d') iv_list(``iv_list_u'') ``opt''
	
	/* Display the result */
	tempname b V N
	matrix `b' = e(b)
	matrix rownames `b' = `varlist'
	matrix colnames `b' = `treatment'
	
	matrix `V' = e(V)
	matrix rownames `V' = `treatment'
	matrix colnames `V' = `treatment'

	scalar `N' = e(N)	
	ereturn post `b' `V', esample(`touse')
	ereturn scalar N = `N'
	ereturn local depvar "`varlist'"
	ereturn local command "``command''"

	ereturn display
end

/* Program to partial out */
capture program drop _partial_out
program define _partial_out
	syntax varlist(min=1 numeric fv ts) [aw] [if] [in], gen(string)
	
	marksample touse
	qui reg `varlist' [`weight'`exp'] if `touse'
	predict `gen', residual, if e(sample)
end


// Unbiased IV estimation with one instrument

/* Main interface */
capture program drop iv_unbiased
program define iv_unbiased, eclass
	syntax varlist(min=1 max=1 numeric) [aw] [if] [in], ///
		treatment(varname) ///
		iv_list(varname) ///
		[vce(string)]
	
	marksample touse
	tempname b b_iv V N XI SIGMA beta flag_failed

	/* Get 1st-stage and reduced-form estimate and variance-covariance matrix */
	_get_xi_sigma `varlist' [`weight'`exp'] if `touse', ///
		treatment(`treatment') iv(`iv_list') vce(`se')

	matrix `b_iv' = e(b_iv)
	matrix `V' = e(V_iv)
	matrix `XI' = e(b)
	matrix `SIGMA' = e(V)

	/* Calculate unbiased IV estimate give xi and sigma */
	_calc_iv_unbiased, xi(`XI') sigma(`SIGMA')
	scalar `flag_failed' = r(flag_failed)
	scalar `beta' = (1 - `flag_failed') * r(beta) + `flag_failed' * `b_iv'[1, 1]
	
	/* Post the estimation result */
	matrix `b' = scalar(`beta')
	matrix rownames `b' = `varlist'
	matrix colnames `b' = `treatment'
	matrix rownames `V' = `treatment'
	matrix colnames `V' = `treatment'

	qui count if `touse'
	scalar `N' = r(N)

	ereturn post `b' `V', esample(`touse')
	ereturn scalar N = `N'
	ereturn scalar flag_failed = `flag_failed'
	ereturn local depvar "`varlist'"
	ereturn local command "iv_unbiased" 
	ereturn display
end

capture program drop _get_xi_sigma
program define _get_xi_sigma, eclass
	syntax varlist(min=1 max=1 numeric) [aw] [if] [in], ///
		treatment(varname) ///
		iv(varname) ///
		[vce(string)]

	marksample touse
	tempname b_iv V_iv XI SIGMA se
	if inlist("`vce'", "", "robust") {
		local `se' "`vce'"
	}
	else if regexm("`vce'", "cluster") {
		local `se' "cluster(`: word 2 of `vce'')"
	}
	else {
		display as error "Unexpected vce specification"
		exit 198
	}
	
	qui ivreg2 `varlist' (`treatment'=`iv') [`weight'`exp'] if `touse', ///
		nocons sfirst savesfprefix(_ivreg2_) ``se''
	
	matrix `b_iv' = e(b)
	matrix `V_iv' = e(V)

	qui est restore _ivreg2_sfirst_`varlist'

	matrix `XI' = e(b)
	matrix `SIGMA' = e(V)
	
	ereturn post `XI' `SIGMA', esample(`touse')
	ereturn matrix b_iv = `b_iv'
	ereturn matrix V_iv = `V_iv'
end


capture program drop _calc_iv_unbiased
program define _calc_iv_unbiased, rclass
	syntax , xi(string) sigma(string)
	
	tempname s2 xi2 tau delta beta flag_failed
	scalar `s2' = sqrt(`sigma'[2, 2])
	scalar `xi2' = `xi'[1, 2] / `s2'
	scalar `tau' = normal(-`xi2') / (normalden(`xi2') * `s2')
	scalar `flag_failed' = missing(`tau')
	if `flag_failed' {
		display as error "WARNING: Cannot calculate inverse Mills ratio", ///
			_newline "Conventional IV estimate is expected to work well"
		scalar `beta' = 0
	}
	else {
		scalar `delta' = `xi'[1, 1] - `xi'[1, 2] * (`sigma'[1, 2] / `sigma'[2, 2])
		scalar `beta' = `tau' * `delta' + (`sigma'[1, 2] / `sigma'[2, 2])
	}
	
	return scalar beta = `beta'
	return scalar flag_failed = `flag_failed'
end


// Unbiased IV estimation with multiple instruments

/* Main interface */
capture program drop iv_unbiased_gmm
program define iv_unbiased_gmm, eclass
	syntax varlist(min=1 max=1 numeric) [aw] [if] [in], ///
		treatment(varname) ///
		iv_list(varlist min=1 numeric) ///
		[nsim(integer 1000000)] ///
		[vce(string)]
	
	marksample touse
	
	/* Get weight and covariance matrices of conventional GMM 2SLS estimation */
	_get_weight_matrix `varlist' [`weight'`exp'] if `touse', ///
		treatment(`treatment') iv_list(`iv_list') vce(`vce')

	tempname V_GMM W_GMM
	matrix `V_GMM' = r(V)
	matrix `W_GMM' = r(W)

	/* Get one-by-one reduced form and 1st-stage estimates */
	_get_xi_sigma_gmm2s `varlist' [`weight'`exp'] if `touse', ///
		treatment(`treatment') iv_list(`iv_list') vce(`vce')

	tempname XI SIGMA 
	matrix `XI' = e(b)
	matrix `SIGMA' = e(V)

	/* Calculate RB estimate via simulation */
	tempname zetas niv ZEROS
	local `niv' = wordcount("`iv_list'")
	local `zetas' ""
	forvalues i = 1(1)``niv'' {
		local `zetas' "``zetas'' zeta1`i' zeta2`i'"
	}
	matrix `ZEROS' = J(1, 2 * ``niv'', 0)

	preserve
		clear
		qui set obs `nsim'

		/* Draw random variables for simulation */
		drawnorm ``zetas'', means(`ZEROS') cov(`SIGMA')

		/* Calculate unbiased estimate */
		forvalues i = 1(1)``niv'' {			
			_gen_unbiased_iv , xi(`XI') sigma(`SIGMA') ///
				zeta1(zeta1`i') zeta2(zeta2`i') pos(`i') gen(beta_`i')
		}
		
		/* Calculate weight on ith unbiased IV estimate */
		_gen_weight , xi(`XI') wt_matrix(`W_GMM') zeta2(zeta2) niv(``niv'') gen(weight)

		/* Calculate RB estimate */
		_calc_RB , b(beta) weight(weight) niv(``niv'')
	restore

	/* Return estimation result */
	tempname b N
	matrix `b' = scalar(r(betaRB))
	matrix rownames `b' = `varlist'
	matrix colnames `b' = `treatment'
	matrix rownames `V_GMM' = `treatment'
	matrix colnames `V_GMM' = `treatment'

	qui count if `touse'
	scalar `N' = r(N)

	ereturn post `b' `V_GMM', esample(`touse')
	ereturn scalar N = `N'
	ereturn local depvar "`varlist'"
	ereturn local command "iv_unbiased_gmm"
	ereturn display
end

/* Program to get weight and covariance matrices in conventional GMM 2SLS */
capture program drop _get_weight_matrix
program define _get_weight_matrix, rclass
	syntax varlist(min=1 max=1 numeric) [aw] [if] [in], ///
		treatment(varname) ///
		iv_list(varlist min=1 numeric) ///
		[vce(string)]
	
	marksample touse
	tempname W V se
	if inlist("`vce'", "", "robust") {
		local `se' "`vce'"
	}
	else if regexm("`vce'", "cluster") {
		local `se' "cluster(`: word 2 of `vce'')"
	}
	else {
		display as error "Unexpected vce specification"
		exit 198
	}
	
	qui ivreg2 `varlist' (`treatment'=`iv_list') [`weight'`exp'] ///
		if `touse', gmm2s noid nocons ``se''
	matrix `W' = e(W)
	matrix `V' = e(V)

	return matrix W = `W'
	return matrix V = `V'
end

/* Program to get one-by-one reduced-form and 1st-stage estimates */
capture program drop _get_xi_sigma_gmm2s
program define _get_xi_sigma_gmm2s, eclass
	syntax varlist(min=1 max=1 numeric) [aw] [if] [in], ///
		treatment(varname) ///
		iv_list(varlist min=1 numeric) ///
		[vce(string)]

	marksample touse
	tempname eq_main instruments i b V

	local `eq_main' ""
	local `instruments' ""
	local `i' = 0
	foreach iv in `iv_list' {
		local `i' = ``i'' + 1
		local `eq_main' "``eq_main'' (eq1``i'': `varlist' - {b1``i''} * `iv')"
		local `eq_main' "``eq_main'' (eq2``i'': `treatment' - {b2``i''} * `iv')"

		local `instruments' "``instruments'' instruments(eq1``i'': `iv', nocons)"
		local `instruments' "``instruments'' instruments(eq2``i'': `iv', nocons)"
	}

	qui gmm ``eq_main'' [`weight'`exp'] if `touse', ///
		``instruments'' `se' winitial(identity)
	matrix `b' = e(b)
	matrix `V' = e(V)
	
	ereturn post `b' `V'
end

/* Program to calculate unbiased IV estimate */
capture program drop _gen_unbiased_iv
program define _gen_unbiased_iv
	syntax , xi(string) sigma(string) zeta1(varname) zeta2(varname) ///
		pos(integer) gen(string)

	tempname s2 s12_s2
	scalar `s2' = sqrt(2 * `sigma'[2*`pos', 2*`pos'])
	scalar `s12_s2' = `sigma'[2*`pos'-1, 2*`pos'] / `sigma'[2*`pos', 2*`pos']
	
	tempvar xi1 xi2 tau delta
	qui gen double `xi1' = `xi'[1, 2*`pos'-1] + `zeta1'
	qui gen double `xi2' = `xi'[1, 2*`pos'] + `zeta2'

	qui gen double `tau' = normal(-`xi2' / `s2') / (normalden(`xi2' / `s2') * `s2')
	qui gen double `delta' = `xi1' - `xi2' * `s12_s2'
	qui gen double `gen' = `tau' * `delta' + `s12_s2'
end

/* Program to calculate weight on unbiased IV estimate */
capture program drop _gen_weight
program define _gen_weight
	syntax , xi(string) wt_matrix(string) zeta2(string) niv(integer) gen(string)

	forvalues i = 1(1)`niv' {
		tempvar xi2`i'
		qui gen double `xi2`i'' = `xi'[1, 2*`i'] - `zeta2'`i'
	}

	tempvar total
	qui gen double `total' = 0
	forvalues i = 1(1)`niv' {
		qui gen double `gen'_`i' = 0
		forvalues j = 1(1)`niv' {
			qui replace `gen'_`i' = `gen'_`i' + `wt_matrix'[`j', `i'] * `xi2`j''
		}
		qui replace `gen'_`i' = `gen'_`i' * `xi2`i''
		qui replace `total' = `total' + `gen'_`i'
	}
	
	forvalues i = 1(1)`niv' {
		qui replace `gen'_`i' = `gen'_`i' / `total'
	}
end

/* Program to calculate RB estimate */
capture program drop _calc_RB
program define _calc_RB, rclass
	syntax , b(string) weight(string) niv(integer)
	
	tempvar b_rb
	qui gen double `b_rb' = 0
		forvalues i = 1(1)`niv' {
			qui replace `b_rb' = `b_rb' + `b'_`i' * `weight'_`i'
	}
	sum `b_rb', meanonly		
	
	tempname betaRB
	scalar `betaRB' = r(mean)
	
	return scalar betaRB = `betaRB'
end
