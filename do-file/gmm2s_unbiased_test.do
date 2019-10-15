clear all
set more off
set matsize 11000

local pwd "D:/GitHub/gmm2s_unbiased"
qui do "`pwd'/do-file/gmm2s_unbiased.do"

/* Program to generate data */
program define _gen_data
	syntax , seed(integer) nobs(integer) mean(string) cov(string) ///
		coef1st(string) coef2nd(string) [clear]
	
	`clear'
	qui set seed `seed'
	qui set obs `nobs'
	
	qui drawnorm x z1 z2 u v, means(`mean') cov(`cov')
	qui gen d = `coef1st'[1, 1] + `coef1st'[1, 2] * x ///
		+ `coef1st'[1, 3] * z1 + `coef1st'[1, 4] * z2 + u
	qui gen y = `coef2nd'[1, 1] + `coef2nd'[1, 2] * x + `coef2nd'[1, 3] * d + v
end

/* Program for progress bar */
capture program drop _update_progress_bar
program define _update_progress_bar
	syntax , iter(integer) [max(integer 50)]
	
	if inrange(`iter', 0.999, 1.001) {
		_set_progress_bar, max(`max')
	}

	if mod(`iter' - 1, `max') == 0 {
		di substr(" " * 4, 1, 4 - strlen(strofreal(`iter' - 1))) _continue

		di strofreal(`iter' - 1) "|" _continue
	}
	di "." _continue

	if mod(`iter', `max') == 0 {
		di _newline
	}
end

capture program drop _set_progress_bar
program define _set_progress_bar
	syntax , [max(integer 50)]
	
	local spaces = " " * 4
	di "`spaces' " _continue
	forvalues j = 1(1)`max' {
		if mod(`j' + 1, 10) == 0 {
			di `j' + 1 _continue
		}
		else if !inrange(mod(`j' + 1, 10), 1, strlen("`=`j'+1'") - 1) { 
			di " " _continue
		}
	}
	di _newline

	di "`spaces'|" _continue
	forvalues j = 1(1)`max' {
		if mod(`j', 10) == 0 {
			di "+" _continue
		}
		else {
			di "-" _continue
		}
	}
	di _newline
end

/* Seed */
local seed = 1120002

/* Number of observations */
local nobs = 1000

/* Number of simulations */
local nsim = 10000000

/* Variance parameters */
foreach var in x z1 z2 u v {
	local s_`var'`var' = 1
}

/* Covariance parameters */
foreach s in s_xz1 s_xz2 s_xu s_xv s_z1z2 s_z1u s_z1v s_z2u s_z2v s_uv {
	local `s' = 0
}
local s_z1z2 = 0.3
local s_uv = 0.5

matrix MEAN = J(1, 5, 0)
matrix COV = [ ///
	`s_xx' , `s_xz1' , `s_xz2' , `s_xu' , `s_xv'  \ ///
	`s_xz1', `s_z1z1', `s_z1z2', `s_z1u', `s_z1v' \ ///
	`s_xz2', `s_z1z2', `s_z2z2', `s_z2u', `s_z2v' \ ///
	`s_xu' , `s_z1u' , `s_z2u' , `s_uu' , `s_uv'  \ ///
	`s_xv' , `s_z1v' , `s_z2v' , `s_uv' , `s_vv' ///
]

/* 1st-stage-equation coefficients */
matrix B1 = [0, 0.2, 0.3, 0.2]

/* 2nd-stage-equation coefficients */
matrix B2 = [0, 0.3, 0.5]

/* Generate data */
_gen_data , seed(`seed') nobs(`nobs') mean(MEAN) cov(COV) ///
	coef1st(B1) coef2nd(B2) clear

/* Estimation */
eststo clear
qui eststo: reg y d, robust
qui eststo: ivreg2 y x (d=z1 z2), gmm2s robust noid
qui eststo: gmm2s_unbiased y, treatment(d) iv_list(z1 z2) controls(x) ///
	nsim(`nsim') vce(robust)

esttab, se nogap label obslast b(%05.4f) se(%05.4f) nonotes nostar ///
	keep(d) mtitle("OLS" "BiasedIV" "UnbiasedIV")


/* Simulation: Number of observations */
local nobs_min = 200
local nobs_max = 20000
local step = 200

local nsim = 100000

forvalues n = `nobs_min'(`step')`nobs_max' {
	_gen_data , seed(`n') nobs(`n') mean(MEAN) cov(COV) ///
		coef1st(B1) coef2nd(B2) clear

	qui gmm2s_unbiased y, treatment(d) iv_list(z1 z2) controls(x) ///
		nsim(`nsim') vce(robust)
	matrix B_unbiased = nullmat(B_unbiased) \ [_b[d]]
	matrix BU_unbiased = nullmat(BU_unbiased) \ [_b[d] + 1.96 * _se[d]]
	matrix BL_unbiased = nullmat(BL_unbiased) \ [_b[d] - 1.96 * _se[d]]

	qui ivreg2 y x (d=z1 z2), gmm2s robust noid
	matrix B_biased = nullmat(B_biased) \ [_b[d]]
	matrix BU_biased = nullmat(BU_biased) \ [_b[d] + 1.96 * _se[d]]
	matrix BL_biased = nullmat(BL_biased) \ [_b[d] - 1.96 * _se[d]]
	
	local iter = `n' / `step'
	_update_progress_bar, iter(`iter')
}

preserve
	foreach stat in B BU BL {
		foreach tag in _unbiased _biased {
			svmat `stat'`tag'
		}
	}
	gen nobs = `step' * _n / 1000
	
	keep if !missing(B_unbiased1)
	twoway (connect B_unbiased1 nobs, color(sky) lp(l)) ///
		(line BU_unbiased1 nobs, lc(sky) lp("-##")) ///
		(line BL_unbiased1 nobs, lc(sky) lp("-##")) ///
		(connect B_biased1 nobs, color(reddish) lp(l)) ///
		(line BU_biased1 nobs, lc(reddish) lp("-##")) ///
		(line BL_biased1 nobs, lc(reddish) lp("-##")), ///
		ytitle("Estimate") ylabel(0(0.1)1, format(%02.1f)) ///
		xtitle("Number of observations (1,000)") xlabel(, nogrid format(%03.1f)) ///
		legend(order(1 "Unbiased IV" 4 "Biased IV") ring(0) pos(2) ///
			region(lw(*0.3) lc(gs13))) ///
		scheme(tt_color)
	graph export "`pwd'/figure/sim_nobs.pdf", replace as(pdf)
		
	save "`pwd'/data/sim_nobs.dta", replace
restore

/* Simulation: Distribution of estimate */
local nrep = 1000
local nobs = 500
local nsim = 100000
capture matrix drop B_unbiased B_biased

forvalues i = 1(1)`nrep' {
	_gen_data , seed(`i') nobs(`nobs') mean(MEAN) cov(COV) ///
		coef1st(B1) coef2nd(B2) clear

	qui gmm2s_unbiased y, treatment(d) iv_list(z1 z2) controls(x) ///
		nsim(`nsim') vce(robust)
	matrix B_unbiased = nullmat(B_unbiased) \ [_b[d]]

	qui ivreg2 y x (d=z1 z2), gmm2s robust noid
	matrix B_biased = nullmat(B_biased) \ [_b[d]]
	
	_update_progress_bar , iter(`i') max(100)
}

preserve
	foreach mat in B_unbiased B_biased {
		svmat `mat'
	}
	
	keep if !missing(B_unbiased1)
	twoway (histogram B_unbiased1, fraction color(sky%60) start(0) width(0.05)) ///
		(histogram B_biased1, fraction color(reddish%60) start(0) width(0.05)), ///
		ytitle("Fraction") ylabel(, format(%03.2f)) ///
		xtitle("Estimate") xlabel(0(0.1)1, nogrid format(%03.2f)) ///
		legend(order(1 "Unbiased IV" 2 "Biased IV") ring(0) pos(2) ///
			region(lw(*0.3) lc(gs13))) ///
		scheme(tt_color)

	graph export "`pwd'/figure/sim_dist.pdf", replace as(pdf)
		
	save "`pwd'/data/sim_dist.dta", replace
restore

/* Simulation: Distribution of estimate (Small sample) */
local nrep = 1000
local nobs = 50
local nsim = 100000
capture matrix drop B_unbiased B_biased

forvalues i = 1(1)`nrep' {
	_gen_data , seed(`i') nobs(`nobs') mean(MEAN) cov(COV) ///
		coef1st(B1) coef2nd(B2) clear

	qui gmm2s_unbiased y, treatment(d) iv_list(z1 z2) controls(x) ///
		nsim(`nsim') vce(robust)
	matrix B_unbiased = nullmat(B_unbiased) \ [_b[d]]

	qui ivreg2 y x (d=z1 z2), gmm2s robust noid
	matrix B_biased = nullmat(B_biased) \ [_b[d]]
	
	_update_progress_bar , iter(`i') max(100)
}

preserve
	foreach mat in B_unbiased B_biased {
		svmat `mat'
	}
	
	keep if !missing(B_unbiased1)
	
	tabstat B_unbiased1 B_biased1, stat(mean sd p1 p10 p25 p50 p75 p90 p99) ///
		columns(stat) format(%04.3f)
	corr B_biased1 B_unbiased1 
	
	twoway (histogram B_unbiased1 if inrange(B_unbiased1, -1, 2), ///
		fraction color(sky%60) start(-1) width(0.05)) ///
		(histogram B_biased1 if inrange(B_biased1, -1, 2), ///
		fraction color(reddish%60) start(-1) width(0.05)), ///
		ytitle("Fraction") ylabel(, format(%03.2f)) ///
		xtitle("Estimate") xlabel(-1(0.25)2, nogrid format(%03.2f)) ///
		legend(order(1 "Unbiased IV" 2 "Biased IV") ring(0) pos(2) ///
			region(lw(*0.3) lc(gs13))) ///
		scheme(tt_color)

	graph export "`pwd'/figure/sim_dist_small.pdf", replace as(pdf)
		
	save "`pwd'/data/sim_dist_small.dta", replace
restore
