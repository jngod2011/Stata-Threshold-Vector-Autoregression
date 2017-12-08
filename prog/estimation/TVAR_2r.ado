// Implement if and in
capture program drop TVAR_2r
program define TVAR_2r, eclass
  syntax varlist(ts min=1) [if] [in], indicator(name) nlag(integer) d(integer) [constant(integer 1) ols(integer 1) no_regime(integer 1)]
  
  //////////////////////////
  if (`nlag' < 1) {
    display as error "Maximum lag order (nlag) inappropropriate"
	exit
  }
  if (`d' > `nlag') {
    display as error "Delay parameter (nlag) inappropropriate"
	exit
  }
  //////////////////////////
  
  ///////////////
  tempname touse
  mark `touse' `if' `in'
  ///////////////
  
  preserve
  keep if `touse'
  display as result ""
  display as result "Threshold Vector Autoregression with 2 Regimes ( TVAR(2) ): Estimation"
  display as text   "Dependent Variables       : `varlist'"
  display as text   "Indicator Variable        : `indicator'"
  display as text   "Delay Parameter (d)       : `d'"
  display as text   "Number of Lags  (nlag)    : 1/`nlag'"
  display as text   "Constant?                 :", _continue
  if (`constant' == 1) {
    display as text "Yes"
  }
  else {
    display as text "No"
  }
  display as text   "Estimation Method         :",_continue
  if (`ols' == 1) {
    display as text "Ordinary Least Square (OLS) Estimation"
  }
  else {
    display as text "Maximum Likelihood (ML) Estimation"
  }
  
  qui TVAR_2r_estimation `varlist', indicator(`indicator') nlag(`nlag') d(`d') constant(`constant') ols(`ols') no_regime(`no_regime')
  
  display as result ""
  
  ////////////////////////////////
  // Header 
  display as result "Estimation Results:"
  display ""
  display as text "{hline 75}"
  display as text "  ", _continue
  display as text %5s  "Regime", _continue
  display as text %10s "Parms", _continue   
  display as text %10s "R^2", _continue     
  display as text %10s "RMSE", _continue
  display as text %15s "Sample Size", _continue
  display as text %15s "Success?"
  display as text "{hline 75}"
  ///////////////////////////////
  display as text %9s  "0", _continue
  local len_varlist = 0
  foreach term in `varlist' {
    local len_varlist = `len_varlist' + 1
  }
  local parm = ((`len_varlist' * `nlag') + `constant')*`len_varlist'
  display as text %10.0f `parm', _continue
  display as text %10.4f e(r0_R_sq), _continue
  display as text %10.4f e(r0_RMSE), _continue
  display as text %15.0f e(r0_size), _continue
  if (e(r0_is_success) == 1) {
    display as text %15s   "Yes"
  }
  else {
    display as text %15s   "No"
  }
  ///////////////////////////////
  display as text %9s "1", _continue
  display as text %10.0f `parm', _continue
  display as text %10.4f e(r1_R_sq), _continue
  display as text %10.4f e(r1_RMSE), _continue
  display as text %15.0f e(r1_size), _continue
  if (e(r1_is_success) == 1) {
    display as text %15s   "Yes"
  }
  else {
    display as text %15s   "No"
  }
  ///////////////////////////////////
  if (`no_regime' == 1) {
    display as text %9s "No Regime", _continue
    display as text %10.0f `parm', _continue
    display as text %10.4f e(no_r_R_sq), _continue
    display as text %10.4f e(no_r_RMSE), _continue
	local s_s = e(r0_size) + e(r1_size)
	display as text %15.0f `s_s', _continue
    if (e(no_r	_is_success) == 1) {
      display as text %15s   "Yes"
    }
    else {
      display as text %15s   "No"
    }
  }
  display as text "{hline 75}"
  //////////////////////////////////
  //////////////////////////////////
  // Header
  display ""
  display as text   "{hline 15}{c TT}{hline 50}"
  display as result "Regime         {c |}  Coefficients "
  display as text   "{hline 15}{c +}{hline 50}"
  /////////////////////////////////////////
  ///////////////////////////////////////////////
  display as text   "Regime 0       {c |}"
  if (e(r0_is_success) == 1) {
  display _column(16) "{c |}"
  foreach dependent in `varlist' {
    display as result %-15s "`dependent':" "{c | }"
	if (`constant' == 1) {
	  local col_num = colnumb( e(r0_beta), "`dependent'_r1:_cons" )
	  display as text   %-15s  "_cons" "{c | }", _continue
	  matrix A = e(r0_beta)
	  local num = A[1,`col_num']
	  display as text   %7.4f `num'
	}
	foreach regressor in `varlist' {
	  forvalues lag = 1/`nlag' {
	    local independent "L`lag'.`regressor'"
		local col_num = colnumb( e(r0_beta), "`dependent'_r1:`independent'_r1" )
	    display as text %-15s "`independent'" "{c | }", _continue
        matrix A = e(r0_beta)
	    local num = A[1,`col_num']
	    display as text  %7.4f `num'
	  }
	}
	display as result _column(16) "{c |}"
  }
  }
  else {
    display as text _column(16)"{c | }Estimation for Regime 0 Unsuccessful"
  }
  display as text   "{hline 15}{c +}{hline 50}"
  //////////////////////////////////////////////////////////////////
  display as text   "Regime 1       {c |}"
  if (e(r1_is_success) == 1) {
  display _column(16) "{c |}"
  foreach dependent in `varlist' {
    display as result %-15s "`dependent':" "{c | }"
	if (`constant' == 1) {
	  local col_num = colnumb( e(r1_beta), "`dependent'_r2:_cons" )
	  display as text   %-15s  "_cons" "{c | }", _continue
	  matrix A = e(r1_beta)
	  local num = A[1,`col_num']
	  display as text   %7.4f `num'
	}
	foreach regressor in `varlist' {
	  forvalues lag = 1/`nlag' {
	    local independent "L`lag'.`regressor'"
		local col_num = colnumb( e(r1_beta), "`dependent'_r2:`independent'_r2" )
	    display as text %-15s "`independent'" "{c | }", _continue
        matrix A = e(r1_beta)
	    local num = A[1,`col_num']
	    display as text  %7.4f `num'
	  }
	}
	display as result _column(16) "{c |}"
  }
  }
  else {
    display as text _column(16)"{c | }Estimation for Regime 1 Unsuccessful"
  }
  if (`no_regime' == 1) {
  display as text   "{hline 15}{c +}{hline 50}"
  //////////////////////////////////////////////////////////////////
  display as text   "No Regime      {c |}"
  if (e(no_r_is_success) == 1) {
  display _column(16) "{c |}"
  foreach dependent in `varlist' {
    display as result %-15s "`dependent':" "{c | }"
	if (`constant' == 1) {
	  local col_num = colnumb( e(no_r_beta), "`dependent':_cons" )
	  display as text   %-15s  "_cons" "{c | }", _continue
	  matrix A = e(no_r_beta)
	  local num = A[1,`col_num']
	  display as text   %7.4f `num'
	}
	foreach regressor in `varlist' {
	  forvalues lag = 1/`nlag' {
	    local independent "L`lag'.`regressor'"
		local col_num = colnumb( e(no_r_beta), "`dependent':`independent'" )
	    display as text %-15s "`independent'" "{c | }", _continue
        matrix A = e(no_r_beta)
	    local num = A[1,`col_num']
	    display as text  %7.4f `num'
	  }
	}
	display as result _column(16) "{c |}"
  }
  }
  else {
    display as text _column(16)"{c | }Estimation for No Regime Unsuccessful"
  }
  }
  display as text   "{hline 15}{c BT}{hline 50}"
  restore
end
