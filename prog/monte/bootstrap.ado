capture program drop han_LR_boot_count_varlist
program define han_LR_boot_count_varlist, rclass
  syntax varlist
  tokenize `varlist'
  while "`*'" != "" {
    local current_var `1'
	local len_varlist = `len_varlist' + 1
    macro shift
  }
  return scalar len = `len_varlist'
end

capture program drop han_LR_boot_check
program define han_LR_boot_check, rclass
  syntax varlist, res(varlist) estimators(name) constant(integer) nsample(integer) initial_obs(integer) [nlag(integer 1)]
  local rt = 0
  dis "in check"
  
  // count: input variables
  han_LR_boot_count_varlist `varlist'
  local len_var = r(len)
  han_LR_boot_count_varlist `res'
  local len_res = r(len)
  //dis "len_res: `len_res'"
  //matrix list `estimators'
  local len_est = colsof(`estimators')
  
  // compute the correct number of res and est based on lags, len_varlist and constant
  if (`len_res' != `len_var') {
    dis as error "Number of residuals does not match with Number of variables"
    local rt = 1
	exit
  }
  
  local num_est = (`constant' + (`len_var' * `nlag'))*`len_var'
  dis "Num estimators: "`num_est'
  dis "Len estimators: "`len_est'
  if (`num_est' != `len_est') {
    display as error "Number of Estimators does not match"
	local rt = 1
	exit
  }
  return scalar is_success = `rt'
  return scalar num_var    = `len_var'
  return scalar num_res    = `len_res'
  return scalar num_est    = `len_est'
  dis "out check"
end

capture program drop han_LR_boot_resampling
program define han_LR_boot_resampling, rclass
  syntax varlist, num_res(integer) [nsample(integer 1)]
  // Re-sampling residuals
  dis "in resampling"
  dis "num_res: `num_res'"
  dis "nsample: `nsample'"

  preserve
	matrix b_res = J(`nsample',`num_res',.)
	local counter = `counter' + 1
	bsample `nsample'
	mkmat `varlist', mat(b_res)
  restore
  //matrix list b_res
  return matrix b_res = b_res
end

capture program drop bootstrap_2r_emp
program define bootstrap_2r_emp, rclass
  //    first  num_var variables : variables
  //    second num_var variables : actual residual
  //    initial_obs means the point of the given fixed point
  syntax varlist, res(varlist) estimators(name) constant(integer) nsample(integer) initial_obs(integer) [nlag(integer 1) set_seed(integer 1234)]
  //local a = colnumb(`estimators',"y1:_cons")
  
  set seed `set_seed' // Set random seed
  
  han_LR_boot_check `varlist', res(`res') estimators(`estimators') constant(`constant') nsample(`nsample') initial_obs(`initial_obs') nlag(`nlag')
  if (r(is_success) == 1) {
    dis "in"
    exit
  }
  local num_var = r(num_var)
  local num_est = r(num_est)
  local num_res = r(num_res)
  
  //dis `num_var'
  //dis `num_est'
  //dis `num_res'
  // Re-sampling residuals
  han_LR_boot_resampling `res', num_res(`num_res') nsample(`nsample')
  matrix b_res = r(b_res)
  matrix list b_res
  
  // Recursively Constructing bootstrapped samples
  preserve
    local counter_dependent = 0
    foreach dependent in `varlist' {
	  local counter_dependent = `counter_dependent' + 1
	  qui gen b_`dependent ' = .
	  // Construct Base Case
	  if (`constant' == 1) {
	    local index_cons = colnumb(`estimators',"`dependent':_cons")
	    replace b_`dependent' = `estimators'[1,`index_cons'] in 1/`nsample'
	  }
	  else {
	    replace b_`dependent' = 0 in 1/`nsample' in 1/`nsample'
	  }
	  // store index of each lag in a matrix
	  // exchange memeory for efficiency
	  matrix index_lags_`dependent' = J(`nlag',`num_var',.)
	  local col_counter = 0
	  
	  foreach regressor in `varlist'  {
	    local col_counter = `col_counter' + 1
	    forvalues nlags = 1/`nlag' {
	      matrix index_lags_`dependent'[`nlags',`col_counter'] = colnumb(`estimators',"`dependent':L`nlags'.`regressor'")
		}
	  }
	  matrix colnames index_lags_`dependent' = `varlist'
	  matrix list index_lags_`dependent'
	  
	  // Construct the first example
	  tokenize `varlist'
	  forvalues nvar = 1/`num_var' {
	    forvalues lag = 1/`nlag' {
		  local estr = `estimators'[1,index_lags_`dependent'[`lag',`nvar']]
		  //dis "estr: `estr' | reg: "`1'[`initial_obs'-`lag']
	      qui replace b_`dependent' = b_`dependent'[1] + `estr' * `1'[`initial_obs'-`lag'] in 1
		}
		qui replace b_`dependent' = b_`dependent'[1] + b_res[1,`nvar'] in 1
		macro shift
	  }
	}
	//dis as result "Subsequent Samples"
	// Recursively construct samples
	local b_sample_var ""
	foreach dependent in `varlist' {
	  local b_sample_var `b_sample_var' b_`dependent'
	  forvalues index = 2/`nsample' {
	    tokenize `varlist'
	    forvalues nvar = 1/`num_var' {
	      forvalues lag = 1/`nlag' {
		    local estr = `estimators'[1,index_lags_`dependent'[`lag',`nvar']]
		    //dis "estr: `estr' | reg: "b_`1'[`initial_obs'-`lag'] 
	        qui replace b_`dependent' = b_`dependent'[`index'] + `estr' * b_`1'[`initial_obs'-`lag'] in `index'
	  	  }
		  qui replace b_`dependent' = b_`dependent'[`index'] + b_res[`index',`nvar'] in `index'
		  macro shift
	    }
	  }
	}
	drop if _n > `nsample'
	//dis "`b_sample_var'"
	mkmat `b_sample_var', matrix(b_sample)
	//matrix colnames b_sample = `varlist'
  restore
  return matrix b_res    = b_res
  return matrix b_sample = b_sample
end
