// TVAR Grid Search

capture program drop TVAR_2r_grid_search
program define TVAR_2r_grid_search, eclass
  syntax varlist(ts min=1) [if] [in], threshold(varname) ptrim(real) d(integer) nlag(integer) [constant(integer 1) ols(integer 1) criterion(integer 1) est_table(integer 0)]
  
  tempname touse
  mark `touse' `if' `in'
  preserve
  keep if `touse'
  
  display ""
  display as result "Threshold Vector Autoregression with 2 Regimes: Grid Search"
  display ""
  
  display as text   "Variables            : `varlist'"
  display as text   "Threshold Variable   : `threshold'"
  display as text   "Delay of Threshold   : "`d'
  display as text   "Number of Lags       : "`nlag'
  display as text   "Estimation Method    :",_continue
  if (`ols' == 1) {
    display as text "Ordinary Least Square Estimation (OLS)"
  }
  else {
    display as text "Maximum Likelihood Estimation (MLE)"
  }
  if (`criterion' == 1) {
    display as text "Grid Search Criterion: Minimum Residual Sum of Square (RSS)"
  }
  else {
    display as text "Grid Search Criterion: Minimum Determinant of the Sigma Matrix"
  }
  dis ""
  
  if (`est_table' == 1) {
    TVAR_2r_grid_op `varlist', threshold(`threshold') ptrim(`ptrim') d(`d') nlag(`nlag') constant(`constant') ols(`ols') criterion(`criterion') est_table(1)
  }
  else {
    qui TVAR_2r_grid_op `varlist', threshold(`threshold') ptrim(`ptrim') d(`d') nlag(`nlag') constant(`constant') ols(`ols') criterion(`criterion') 
  }
  //ereturn list
  
  dis ""
  display as result "Grid Search Results:"
  if (e(num_success) >= 1) {
  display as text   "{hline 15}{hline 55}"
  if (`criterion' == 1) {
    display as text   "     Min RSS      Threshold     Index           Num of "
	display as text   "                                           Successful Trials"
	display as text   "{hline 15}{hline 55}"
  
    display as result "  " %10.0g e(minStat), _continue
    display as result "  " %12.0g e(minStat_thrd), _continue
    display as result "  " %7.0g e(minStat_thrd_index), _continue
    display as result "  " %20.0g e(num_success)
  }
  else {
    display as text   "     Min Det      Threshold     Index           Num of "
	display as text   "                                           Successful Trials"
	display as text   "{hline 15}{hline 55}"
	
	display as result "  " %10.0g e(minStat), _continue
    display as result "  " %12.0g e(minStat_thrd), _continue
    display as result "  " %7.0g e(minStat_thrd_index), _continue
    display as result "  " %20.0g e(num_success)
  }
  //matrix list e(minDet_r1_beta)
  dis ""
  local column ""
  
  local len_varlist = 0
  foreach dependent in `varlist' {
    local len_varlist = `len_varlist' + 1
    local column `column' `dependent'
    foreach regressor in `varlist' {
	  forvalues lag = 1/`nlag' {
	    local column `column' L`lag'.`regressor'
	  }
	}
	if (`constant' == 1) {
	  local column `column' "_cons"
	}
  }
  }
  else {
    display as text "No trials are successful."
  }
  restore
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
end
