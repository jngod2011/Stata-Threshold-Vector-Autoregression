capture mata: mata drop TVAR_2r_estimation_cons_tol_res() 
capture mata: mata drop TVAR_est_cons_Sigma()

capture program drop TVAR_2r_estimation
program define TVAR_2r_estimation, eclass
  // last variable being the indicator
  syntax varlist(ts min=1), indicator(name) nlag(integer) d(integer) [constant(integer 1) ols(integer 1) no_regime(integer 0)]

  // count the number of variables 
  local vars `varlist'
  tokenize `vars'
  local sample_size = _N
  local length_varlist = 0
  while "`*'" != "" {
	local length_varlist = `length_varlist' + 1
    macro shift
  }

  // Check if d appropriate
  if (`d' > `nlag') {
    display as error "Delay Parameter (d) must be smaller than or equal to the lag order."
    exit
  }
  
  // Check if indicator valid
  local sample_size = _N
  local indicator_size = `sample_size' - `d'
  
  preserve
    // Classify the observations into 2 subsets by indicators
    local vars `varlist'
    tokenize `vars'
	local counter = 0
	local r1_varlist_for_var ""
	local r2_varlist_for_var ""
    while "`*'" != "" {
	  local current_var `1'
	  qui gen `current_var'_r1 = .
	  qui gen `current_var'_r2 = .
	  if ("`current_var'" != "`indicator'") {
	    local r1_varlist_for_var `r1_varlist_for_var' `current_var'_r1
	    local r2_varlist_for_var `r2_varlist_for_var' `current_var'_r2
	  }
	  local r1_counter = 0
	  local r2_counter = 0
	  
	  forvalues i = 1/`indicator_size' {
	    if (`indicator'[`i'] == 0) {
		  // Regime 1
		  local r1_counter = `r1_counter' + 1
		  qui replace `current_var'_r1 = `current_var'[`i'] in `r1_counter'
		}
		else {
		  // Regime 2
		  local r2_counter = `r2_counter' + 1
		  qui replace `current_var'_r2 = `current_var'[`i'] in `r2_counter'
		}
	  }
	  macro shift
	}
	//display "R1: "`r1_counter'
	//display "R2: "`r2_counter'
	list `r1_varlist_for_var'
	// Run VAR on r1
	display as text ""
    if (`constant' == 0) {
	  if (`ols' == 1) { 
	    display as text "OLS  Estimation: Regime 0  (No Constant)"
	    var_ols `r1_varlist_for_var', constant(0) nlag(`nlag') sample_size(`r1_counter')
	  }
	  else {
	    display as text "ML  Estimation: Regime 0  (No Constant)"
	    var_mle `r1_varlist_for_var', constant(0) nlag(`nlag') sample_size(`r1_counter')
	  }
	}
	else {
	  if (`ols' == 1) {
	    display as text "OLS  Estimation: Regime 0  (With Constant)"
	    var_ols `r1_varlist_for_var', nlag(`nlag') sample_size(`r1_counter')
	  }
	  else {
	    display as text "ML  Estimation: Regime 0  (With Constant)"
	    var_mle `r1_varlist_for_var', nlag(`nlag') sample_size(`r1_counter')
	  }
	}
	// Collect statistics from running VAR on r1
	local r0_is_success = e(is_success)
	local  r0_df         = e(df)
	if (`r0_is_success' == 1) {
	  local  r0_trace      = e(Sigma_trace)
	  local  r0_det        = e(Sigma_det)
	  local  r0_RMSE       = e(RMSE)
	  local  r0_R_sq       = e(R_sq)
	  matrix r0_Sigma    = e(Sigma)
	  matrix r0_res        = e(res)
	  matrix r0_beta       = e(beta)
	}
	////////////// End of r1 ///////////////////
	
	// Run VAR on Regime 2
	display ""
	if (`constant' == 0) {
	  if (`ols' == 1) {
	    display as text "OLS Estimation: Regime 1  (No Constant)"
	    var_ols `r2_varlist_for_var', constant(0) nlag(`nlag') sample_size(`r2_counter')
	  }
	  else {
	    display as text "ML  Estimation: Regime 1  (No Constant)"
	    var_mle `r2_varlist_for_var', constant(0) nlag(`nlag') sample_size(`r2_counter')
	  }
	}
	else {
	  if (`ols' == 1) {
	    display as text "OLS  Estimation: Regime 1  (With Constant)"
	    var_ols `r2_varlist_for_var', nlag(`nlag') sample_size(`r2_counter')
	  }
	  else {
	    display as text "ML  Estimation: Regime 1  (With Constant)"
	    var_mle `r2_varlist_for_var', nlag(`nlag') sample_size(`r2_counter')
	  }
	}
	// Collect Statistics for r1
	local r1_is_success = e(is_success)
	local  r1_df         = e(df)
	if (`r1_is_success' == 1) {
	  local  r1_trace      = e(Sigma_trace)
	  local  r1_det        = e(Sigma_det)
	  local  r1_RMSE       = e(RMSE)
	  local  r1_R_sq       = e(R_sq)
	  matrix r1_Sigma    = e(Sigma)
	  matrix r1_res        = e(res)
	  matrix r1_beta       = e(beta)
	}
	
	// No regime
	if (`no_regime' == 1) {
	  display ""
	  if (`constant' == 0) {
	    if (`ols' == 1) {
	      display as text "OLS  Estimation: No Regime (No Constant)"
	      var_ols `varlist', constant(0) nlag(`nlag')
	    }
	    else {
	      display as text "ML  Estimation: No Regime (No Constant)"
	      var_mle `varlist', constant(0) nlag(`nlag')
	    }
	  }
	  else {
	    if (`ols' == 1) {
	      display as text "OLS  Estimation: No Regime (With Constant)"
	      var_ols `varlist', nlag(`nlag')
	    }
	    else {
	      display as text "ML  Estimation: No Regime (With Constant)"
	      var_mle `varlist', nlag(`nlag')
	    }
	  }
	  // Collect statistics for no regime
	  local no_r_is_success = e(is_success)
	  local  no_r_df         = e(df)
	  if (`no_r_is_success' == 1) {
	    local  no_r_trace      = e(Sigma_trace)
	    local  no_r_det        = e(Sigma_det)
	    local  no_r_RMSE       = e(RMSE)
	    local  no_r_R_sq       = e(R_sq)
	    matrix no_r_Sigma      = e(Sigma)
	    matrix no_r_res        = e(res)
	    matrix no_r_beta       = e(beta)
	  }
	}
	// Construct total residual (for r1 and r2)
	qui putmata indicator = `indicator', replace
	mata: r0_res = st_matrix("r0_res");
	if (`r0_is_success' != 1) {
	  mata: r0_res = J(`r0_df',1,.)
	}
	mata: r1_res = st_matrix("r1_res");
	if (`r1_is_success' != 1) {
	  mata: r1_res = J(`r1_df',1,.)
	}
	if (`r1_is_success' == 1 & `r0_is_success' == 1) {
	  mata: TVAR_2r_estimation_tol_res = TVAR_2r_estimation_cons_tol_res(r0_res, r1_res, indicator, `sample_size');
	  //dis "Tol Res:"
	  //mata: TVAR_2r_estimation_tol_res;
	
	  // Construct the total var-cov matrix
	  mata: TVAR_2r_estimation_Sigma = TVAR_est_cons_Sigma (TVAR_2r_estimation_tol_res);
	  mata: st_matrix("r_Sigma",TVAR_2r_estimation_Sigma);
	  mata: st_matrix("r_det",det(TVAR_2r_estimation_Sigma));
	  mata: st_matrix("r_trace",trace(TVAR_2r_estimation_Sigma));
	}
  restore
  ////////////////////// End of Estimation //////////////////////////
  ereturn clear
  ereturn scalar r0_is_success = `r0_is_success'
  ereturn scalar r0_df         = `r0_df'
  ereturn scalar r0_size       = `r1_counter'
  if (`r0_is_success' == 1) {
    ereturn scalar r0_trace      = `r0_trace'
	ereturn scalar r0_det        = `r0_det'
	ereturn scalar r0_RMSE       = `r0_RMSE'
	ereturn scalar r0_R_sq       = `r0_R_sq'
	ereturn matrix r0_Sigma    = r0_Sigma
	ereturn matrix r0_res        = r0_res
	ereturn matrix r0_beta       = r0_beta
  }

  ereturn scalar r1_is_success = `r1_is_success'
  ereturn scalar r1_df         = `r1_df'
  ereturn scalar r1_size       = `r2_counter'
  if (`r1_is_success' == 1) {
    ereturn scalar r1_trace      = `r1_trace'
	ereturn scalar r1_det        = `r1_det'
	ereturn scalar r1_RMSE       = `r1_RMSE'
	ereturn scalar r1_R_sq       = `r1_R_sq'
	ereturn matrix r1_Sigma    = r1_Sigma
	ereturn matrix r1_res        = r1_res
	ereturn matrix r1_beta       = r1_beta
  }
  
  if (`r0_is_success' == 1 & `r1_is_success' == 1) {
    ereturn matrix r_Sigma = r_Sigma
	//display "trace: "r_trace[1,1]
	ereturn scalar r_Sigma_det     = r_det[1,1]
	ereturn scalar r_Sigma_trace   = r_trace[1,1]
  }
  
  if (`no_regime' == 1) {
    ereturn scalar no_r_is_success = `no_r_is_success'
    ereturn scalar no_r_df         = `no_r_df'
    if (`no_r_is_success' == 1) {
      ereturn scalar no_r_trace      = `no_r_trace'
      ereturn scalar no_r_det        = `no_r_det'
	  ereturn scalar no_r_RMSE       = `no_r_RMSE'
	  ereturn scalar no_r_R_sq       = `no_r_R_sq'
	  ereturn matrix no_r_Sigma    = no_r_Sigma
	  ereturn matrix no_r_res        = no_r_res
	  ereturn matrix no_r_beta       = no_r_beta
    }
  }
end

mata:

real matrix TVAR_est_cons_Sigma(tol_res) {
  real matrix Sigma, tmp1, tmp2;
  num_cols = cols(tol_res);
  Sigma = J(num_cols, num_cols, .);
  
  for (a = 1; a <= num_cols; a++) {
    tmp1 = (tol_res)[.,a];
    for (b = 1; b <= num_cols; b++) {
	  tmp2 = (tol_res)[.,b];
	  Sigma[a,b] = (tmp1') * tmp2;
	}
  }
  return (Sigma);
}

real matrix TVAR_2r_estimation_cons_tol_res(r1_res, r2_res, ind, sample_size) {
  real matrix tol_res;
  //printf("A\n");
  num_rows_r1 = rows(r1_res);
  num_rows_r2 = rows(r2_res);
  num_rows = num_rows_r1 + num_rows_r2;
  num_cols = cols(r1_res);
  tol_res = J(num_rows,num_cols,.);
  num_rows_ind = rows(ind);
  
  real scalar r1_counter, r2_counter;
  
  start_ind_index = -(num_rows - sample_size) + 1;
  /*
  printf("Ind: \n");
  ind;
  printf("r1_res: \n");
  r1_res;
  printf("r2_res: \n");
  r2_res;
  
  printf("Row(tol_res)    = %f \n",num_rows);
  printf("Row(r1_res)     = %f \n",num_rows_r1);
  printf("Row(r2_res)     = %f \n",num_rows_r2);
  printf("Row(ind)        = %f \n",num_rows_ind);
  printf("Sample Size     = %f \n",sample_size);
  printf("Start ind index = %f \n",start_ind_index);
  printf("Col             = %f \n",num_cols);
  */
  for (ncol = 1; ncol <= num_cols; ncol++) {
    r1_counter = 0; r2_counter = 0; tol_counter = 0
    for (row = sample_size; row >= 1; row--) { // ind_counter
	  if (ind[row,1] == 0 & r1_counter < num_rows_r1) {
	    r1_counter++;
		tol_counter++;
	    tol_res[tol_counter,ncol] = r1_res[r1_counter,ncol];
	  }
	  if (ind[row,1] == 1 & r2_counter < num_rows_r2) {
	    r2_counter++;
		tol_counter++;
	    //tol_res[tol_counter,ncol] = 4;
		tol_res[tol_counter,ncol] = r2_res[r2_counter,ncol];
	  }
	  /*
	  printf("Row: %f ", row);
	  printf("Col: %f | ", ncol);
	  printf("tol: %f | ", tol_counter);
	  printf("ind: %f | ", ind[row,1]);
	  printf("R1 : %f | ", r1_counter);
	  printf("R2 : %f \n ", r2_counter);
	  */
	  //tol_res[tol_counter,ncol] = ind[row,1];
	  //printf("r1: %f | r2: %f | tol: %f \n",r1_counter,r2_counter,tol_counter);
	}
	//printf("R1: %f | R2: %f \n",r1_counter,r2_counter);
  }
  //printf("TVAR_est_cons_Sigma out of loop\n");
  //tol_res;
  return (tol_res);
}

end
