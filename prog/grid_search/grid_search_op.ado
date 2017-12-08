/************************************************
 Purpose: Run a Threshold Vector Autoregression
          with 2 Regimes with Grid Search
************************************************/

capture mata: mata drop find_min_stat() 
capture mata: mata drop cons_indicator() 
capture mata: mata drop cons_y_d() 
capture mata: mata drop cons_threshold() 
capture mata: mata drop order_confirm()

// Remember to enforce dependency among programmes
capture program drop TVAR_2r_grid_op
program define TVAR_2r_grid_op, eclass
  syntax varlist(ts min=1), threshold(varname) ptrim(real) d(integer) nlag(integer) [constant(integer 1) ols(integer 1) criterion(integer 1) est_table(integer 0)]
  
  local num_r1 = 0
  local num_r2 = 0
  local num_success = 0
  
  // Count how many variables in the varlist
  local vars `varlist'
  tokenize `vars'
  local sample_size = _N
  local length_varlist = 0
  local threshold_exist = 0
  while "`*'" != "" {
	if ("`1'" == "`threshold'") {
	  local threshold_exist = 1
	}
	local length_varlist = `length_varlist' + 1
    macro shift
  }
  
  // Check erratic options
  if (`ptrim' > 1 | `ptrim' < 0) {
    exit
  }
  if (`d' > `sample_size') {
    exit
  }
  if (`nlag' > `sample_size') {
    exit
  }
  if (`d' > `nlag') {
    exit
  }
  
  // Generate threshold
    qui putmata dependent = (`varlist'), replace
	qui putmata threshold = (`threshold'), replace
    // 1. Create the threshold vector
	mata: _gamma   = cons_threshold(threshold,`ptrim');
	mata: st_matrix("_gamma",_gamma);
	mata: size_gamma = rows(_gamma);
	mata: y_d      = cons_y_d(threshold,`d');
	mata: size_y_d = rows(y_d);
	mata: st_matrix("size_y_d",size_y_d);
	mata: st_matrix("size_gamma_",size_gamma);
	
  // Generate indicator function
  preserve
	// 2. Generate y_d (threshold determinant)
	local size_yd = size_y_d[1,1]

	// 2. Estimation
	local size_gamma = size_gamma_[1,1]
	
	// 1st column: threshold
	// 2nd column: threshold index
	// 3rd column: RSS
	// 4th column: Det
	// 5th column: is_success
	matrix TVAR_stat = J(`size_gamma',5,.)
	matrix colnames TVAR_stat = Threshold Thrd_Index RSS Det is_success
	local gamma_index = 1
	// Interface
	display ""
	display as text   "Number of Trials  : "`size_gamma'

	/////////////////////////////
	local gamma_index = 1
	dis ""
	forvalues gamma_index = 1/`size_gamma' {
	//forvalues gamma_index = 1/40 {
	  
	  display as result `gamma_index ', _continue
	  if (mod(`gamma_index',10) == 0) {
	    dis ""
	  }
	  
	  // construct indicator corresponding to gamma
	  local threshold_val = _gamma[`gamma_index',1]
	  mata: TVAR_indicator = cons_indicator(y_d,`threshold_val',`d');
	  qui getmata indicator = TVAR_indicator, replace
	  capture TVAR_2r_estimation `varlist', indicator(indicator) nlag(`nlag') d(`d') constant(`constant') ols(`ols')
	  qui drop indicator
	  //ereturn list
	  
	  matrix TVAR_stat[`gamma_index',1] = `threshold_val'
	  matrix TVAR_stat[`gamma_index',2] = _gamma[`gamma_index',2]
	  if (e(r0_is_success) == 1 & e(r1_is_success) == 1) {
	    local num_success = `num_success' + 1
	    matrix TVAR_stat[`gamma_index',5] = 1
	    matrix TVAR_stat[`gamma_index',3] = e(r_Sigma_trace)
	    matrix TVAR_stat[`gamma_index',4] = e(r_Sigma_det)
	  }
	  else {
	    matrix TVAR_stat[`gamma_index',5] = 0
	  }
	  //matrix list TVAR_stat
    }
	dis ""
	//matrix list TVAR_stat
  if (`num_success' >= 1) {
  ////////////////////////////////////
  // Find min RSS
  //matrix list TVAR_stat
  mata: TVAR_stat = st_matrix("TVAR_stat");
  if (`criterion' == 1) {
    mata: criterion = 3;
	//ereturn string criterion "Minimum RSS"
  }
  else {
    mata: criterion = 4;
	//ereturn scalar criterion "Minimum Det"
  }
  mata: is_valid = 5;
  //mata: printf("Criterion: %f | is_valid: %f\n",criterion,is_valid);
  mata: _min_index = find_min_stat(TVAR_stat, criterion, is_valid);
  mata: st_matrix("_min_index",_min_index);
  //dis "Hello: "_min_index[1,1]
  
  local min_index = TVAR_stat[_min_index[1,1],2]
  if (`criterion' == 1) {
    local minStat       = TVAR_stat[_min_index[1,1],3]
  }
  else {
    local minStat       = TVAR_stat[_min_index[1,1],4]
  }
  local minStat_thrd  = TVAR_stat[_min_index[1,1],1]
  local minStat_thrd_index  = TVAR_stat[_min_index[1,1],2]
  ////////////////////////////////////
  
  ////////////////////////////////////
  // Statistics of the one with min RSS
  // construct indicator corresponding to gamma
  //quietly{
  local threshold_val = `minStat_thrd'
  //dis "threshold_val: "`threshold_val'
  mata: indicator = cons_indicator(y_d,`threshold_val',`d');
  mata: st_matrix("minStat_indicator",indicator);
  qui getmata indicator = indicator, force replace
  
  
  if (`est_table' == 0) {
    TVAR_2r_estimation `varlist', indicator(indicator) nlag(`nlag') d(`d') constant(`constant') ols(`ols') no_regime(1)
  }
  else {
    TVAR_2r `varlist', indicator(indicator) nlag(`nlag') d(`d') constant(`constant') ols(`ols') no_regime(1)
  }
  
  
  qui drop indicator
  //ereturn list
  
  matrix minStat_r0_beta = e(r0_beta)
  //matrix list minRSS_r1_beta
  matrix minStat_r1_beta = e(r1_beta)
  //matrix list minRSS_r2_beta
  matrix minStat_no_r_beta = e(no_r_beta)
  
  matrix minStat_r0_Sigma = e(r0_Sigma)
  matrix minStat_r1_Sigma = e(r1_Sigma)
  matrix minStat_r_Sigma  = e(r_Sigma)
  
  ///////////////////////////////////
  }
  restore
  ////////////////////////////////////
  // Return
  ereturn clear
  ereturn scalar num_success        = `num_success'
  if (`num_success' >= 1) {
  ereturn scalar minStat            = `minStat'
  ereturn scalar minStat_thrd       = `minStat_thrd'
  ereturn scalar minStat_thrd_index = `minStat_thrd_index'
  
  ereturn matrix minStat_r0_beta    = minStat_r0_beta
  ereturn matrix minStat_r1_beta    = minStat_r1_beta
  ereturn matrix minStat_no_r_beta  = minStat_no_r_beta
  ereturn matrix minStat_r0_Sigma  = minStat_r0_Sigma
  ereturn matrix minStat_r1_Sigma  = minStat_r1_Sigma
  ereturn matrix minStat_r_Sigma   = minStat_r_Sigma
  ereturn matrix Stat              = TVAR_stat
  ereturn matrix minStat_indicator = minStat_indicator
  }
end

mata:
mata set matastrict off

real scalar find_min_stat(TVAR_stat, cri, is_valid) {
  min_index = 0;
  num_rows = rows(TVAR_stat);
  //printf("in here\n");
  for (i = 1; i <= num_rows; i++) {
    if (min_index == 0) {
	  if (TVAR_stat[i,is_valid] != 1) {
	    do {
		  i++;
		} while (TVAR_stat[i,is_valid] != 1);
	  }
	  //printf("before\n");
	  min_index = i;
	  //printf("pre-set min_index: %f\n",i);
	}
	else {
	  if (TVAR_stat[i,is_valid] == 1) {
	    //printf("in Here A\n");
	    if (TVAR_stat[i,cri] < TVAR_stat[min_index,cri]) {
		  //printf("in Here\n");
		  min_index = i;
		}
	  }
	}
  }
  return (min_index);
}

real matrix cons_indicator(y_d, _g, _d) {
  real matrix indicator;
  num_rows = rows(y_d);
  indicator = J(num_rows,1,.);
   
  for (t = 1; t <= num_rows; ++t) {
    if (t - _d > 0) {
	  //printf("t: %f | yd: %f | g: %f\n",t,y_d[t,1],_g);
      if (y_d[t,1] <= _g) {
	    indicator[t,1] = 1;
	  }
	  else {
	    indicator[t,1] = 0;
	  }
	}
  }
  return (indicator);
}

real matrix cons_y_d(y, _d) {
  // 2. Generate y_{t - d}
  num_rows = rows(y);
  start_index = (_d+1);
  end_index   = num_rows;
  //printf("d: %f | %f - %f | %f\n",_d,start_index,end_index,rows(y));
  //printf("y_d\n");
  y_d = J(num_rows, 1, .);
  for (t = 1; t <= num_rows; t++) {
    if (t - _d >= 1) {
	  y_d[t,1] = y[t - _d,1];
	}
  }
  //y_d = y[(_d+1)::num_rows,.];
  return (y_d);
}

real matrix cons_threshold(threshold, ptrim) {
  num_rows = rows(threshold);
  ntrim    = round(num_rows * ptrim);
  start_index    = ntrim + 1;
  end_index      = start_index + num_rows - (ntrim * 2) - 1;
  real matrix trimmed_sorted_threshold;
  trimmed_sorted_threshold = threshold;
  //printf("num_rows: %f | ntrim: %f | start: %f | end: %f \n",num_rows,ntrim,start_index,end_index);
  
  // 1. Create a vector of index and merge it with threshold
  index = J(num_rows,1,.);
  for (i = 1; i <= num_rows; i++) {
    index[i,1] = i;
  }
  trimmed_sorted_threshold = threshold,index;
  //threshold = threshold[start_index..end_index,.];
  // Sort threshold
  trimmed_sorted_threshold = sort(trimmed_sorted_threshold,1);
  //printf("Start: %f | End: %f \n",start_index,end_index);
  trimmed_sorted_threshold = trimmed_sorted_threshold[start_index..end_index,.];
  //threshold = order_confirm(threshold, 1);
  return (trimmed_sorted_threshold);
}

real matrix function order_confirm(threshold, criterion) {
   size = rows(threshold);
   if (threshold[1,1] < threshold[2,1]) {
     return(threshold)
   }
   if (threshold[1,1] == threshold[2,1]) {
     // trace till order identified
	 size = rows(threshold);
	 for (i=3; i <= size; ++i) {
	   if (threshold[1,1] < threshold[i,1]) {
	     return(threshold);
	   }
	 }
   }
   // descending - flip the matrix to ascending
   flip_threshold = J(size,2,.);
   for(i = 1; i <= size; ++i) {
     flip_threshold[i,1] = threshold[size-i+1,1];
	 flip_threshold[i,2] = threshold[size-i+1,2];
   }
   return(flip_threshold); 
}
end
