capture mata: mata drop cons_beta() cons_est() cons_critical_val()
capture program drop sim_2r_normal
program define sim_2r_normal, rclass
  args num_regressors b_times b_nsample is_constant nlag d initial_obs ptrim lv_sig
  //syntax num(integer),num_regressors(integer) //b_times(integer) b_nsample(integer) is_constant(integer) nlag(integer) initial_obs(integer)
  if ("`num_regressors'" == "" | "`b_times'" == "" | "`b_nsample'" == "" | /*
  */  "`is_constant'" == "" | "`nlag'" == "" | "`initial_obs'" == "") {
    display as error "Program input missing."
	exit
  }
  set more off
  
  ///////////////////////////////
  // Pre-Table
  display as result "Threshold Vector Autoregression : Monte Carlo Simulation"
  display as text   "Bootstrap Sample Size           : "`b_nsample'
  display as text   "Number of Bootstrapping         : " `b_times'
  display as text   "Initial Index of Observation    : "`initial_obs'
  display ""
  display as text   "Number of Uniform[0,1] Variables: 1 to "`num_regressors'
  display as text   "Number of Normal (0,1) Residuals: 1 to "`num_regressors'
  display as text   "Number of Lags                  : "`is_constant'
  display as text   "Threshold Delay                 : "`d'
  display as text   "Percentage of Trimming          : "`ptrim'
  
  ///////////////////////////////
  matrix det_mat     = J(`b_times',`num_regressors',.)
  local num_coef = (`is_constant' + (`num_regressors' * `nlag'))*`num_regressors'
  mata: beta_boot_normal = cons_beta(`num_coef');
  
  local num_reg = 1
  local set_seed = 1
  //forvalues b = 1/1 {
  forvalues b = 1/`b_times' {
    display ""
    display "Simulation time: "`b'
    display ""
	forvalues num_reg = 1/`num_regressors' {
    //forvalues num_reg = 1/2 {
      display "Number of regressors: ",`num_reg'
      preserve
        local varlist ""
	    local reslist ""
	
      forvalues i = 1/`num_reg' {
	    local varlist `varlist' var_`i'
	    local reslist `reslist' res_`i'
	    qui gen var_`i' = runiform()
	    qui gen res_`i' = rnormal()
	  }
	  //list `varlist'
	  local beta_index = (`is_constant' + (`num_reg' * `nlag'))*`num_reg'
	  //display "Index: "`beta_index'
	  //mata: beta_boot_normal;
	  mata: _est_ = beta_boot_normal[1,1::`beta_index'];
	  //mata: _est_ = cons_est(beta_boot_normal,`beta_index');
	  mata: st_matrix("_est_",_est_);
	  /////////////////////////////////////
	  // name columns
      local est_header ""
      foreach dep in `varlist' {
        if (`is_constant' == 1) {
	      local est_header `est_header' "`dep':_cons"
	    }
        foreach regressor in `varlist' {
	      forvalues i = 1/`nlag' {
	        local est_header `est_header' "`dep':L`i'.`regressor'"
	      }
	    }
      }
      matrix colnames _est_ = `est_header'
	  /////////////////////////////////////
	  // Bootstrapping
	  local set_seed = `b' * `num_reg'
	  qui bootstrap_2r_emp `varlist', res(`reslist') estimators(_est_) constant(`is_constant') nsample(`b_nsample') initial_obs(3) nlag(1) set_seed(`set_seed')
	  ////////////////////////////////////
	  // Estimations
	  TVAR_2r_grid_op `varlist', threshold(var_1) ptrim(`ptrim') d(`d') nlag(`nlag') constant(`is_constant') ols(1)
	
	  matrix minStat_r_Sigma = e(minStat_r_Sigma)
	  
	  local stat = log(det(minStat_r_Sigma))
	  //local stat = log(e(r1_det) + e(r2_det)) // MLE
	  //dis "Stat: "`stat'
	  //dis "
	  if (`is_constant' == 0) {
	    qui var_ols `varlist', nlag(`nlag') constant(0)
	  }
	  else {
	    qui var_ols `varlist', nlag(`nlag') constant(1)
	  }
	  //ereturn list
	  //dis "after mle"
	  matrix A = e(Sigma)
	  local sum_det = det(A)
	  local stat = (log(`sum_det') - `stat') * `b_nsample'
	  matrix det_mat[`b',`num_reg'] = `stat'
	  
	  ///////////////////////////////////
	drop `varlist'
    restore
  }
}
mata: boot_normal_det_mat = st_matrix("det_mat");
mata: cval = cons_critical_val(boot_normal_det_mat,`lv_sig');
mata: st_matrix("cval",cval);
//mata: cval;
local _width = `num_regressors' * 10 + 20

local percentage__ = (1 - `lv_sig') * 100
display as result "`percentage__'% Critical Values:"

display as text ""
display as text "{hline `_width'}"
display as text "#Regressors  ", _continue
// Header
forvalues i = 1/`num_regressors' {
  display as text %10.0f `i', _continue
}
display as text ""
display as text "{hline `_width'}"
display as text "             ", _continue
forvalues i = 1/`num_regressors' {
  display as result %10.4f cval[1,`i'], _continue
}
display as text ""
display as text "{hline `_width'}"

// Return
return matrix test_stat = det_mat
return matrix lv_sig_critical_val = cval
end

mata: 
real matrix cons_beta(size) {
  real matrix beta;
  rseed(1234)
  beta = runiform(1,size);
  return (beta);
}

real matrix cons_est(beta,beta_index) {
  real matrix _est_;
  num_rows = rows(beta);
  //print_matrix(&beta)
  _est_ = J(num_rows,beta_index,.);
  //print_matrix(&_est_)
  for (col = 1; col <= beta_index; col++) {
    _est_[1,col] = beta[1,col];
  }
  //print_matrix(&_est_);
  return (_est_);
}

real matrix cons_critical_val (det_mat,lv_sig) {
  num_col = cols(det_mat);
  num_row = rows(det_mat);
  real matrix colmat, cr_mat;
  //lv_sig = 0.05;
  pos = num_row - round((num_row + 1) * lv_sig) + 1;
  //printf("Pos: %f\n",pos);
  cr_mat = J(1,num_col,.)
  for (col = 1; col <= num_col; col++) {
    colmat = det_mat[.,col];
	colmat = sort(colmat,1);
	cr_mat[1,col] = colmat[pos,1];
  }
  return (cr_mat);
}


end

/*

matrix list det_mat
mata: det_mat      = st_matrix("det_mat");

mata:
real matrix cons_critical_val (det_mat) {
  num_col = cols(det_mat);
  num_row = rows(det_mat);
  real matrix colmat, cr_mat;
  lv_sig = 0.05;
  pos = num_row - round((num_row + 1) * lv_sig) + 1;
  cr_mat = J(1,num_col,.)
  for (col = 1; col <= num_col; col++) {
    colmat = det_mat[.,col];
	colmat = sort(colmat,1);
	cr_mat[1,col] = colmat[pos,1];
  }
  return (cr_mat);
}
end
*/
