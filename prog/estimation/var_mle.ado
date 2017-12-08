capture program drop var_mle
program define var_mle, eclass
  syntax varlist [, constant(integer 1) nlag(integer 1) sample_size(integer -1)]
  
  local is_success = 1
  
  if (`sample_size' == -1) {
    local sample_size = _N
  }
  //dis "Sample size: `sample_size'"
  // Count how many variables in the varlist
  local vars `varlist'
  tokenize `vars'
  local length_varlist = 0
  while "`*'" != "" {
    local current_var `1'
	local length_varlist = `length_varlist' + 1
    macro shift
  }
  
  // Run VAR MLE
  preserve
    qui drop if _n > `sample_size'
    if (`constant' == 1) {
      qui capture var `varlist', lags(1/`nlag')
    }
    else {
      qui capture var `varlist', lags(1/`nlag') noconstant
    }
	if (!_rc) {
	  // list `varlist'
	  qui capture predict res, res
	  if (_rc) {
	    display as error "Error in MLE Prediction"
		local is_success = 0
		ereturn scalar is_success = 0
		exit
	  }
	  
	  local start_ind = 1 + `nlag'
	  // Export residual matrices
	  local _size_ = `sample_size' - `start_ind' + 1
	  local is_first = 0
	  
	  foreach term in `varlist' {
	    // res
		qui predict `term'_r, res equation(`term') 
		qui putmata `term'_r = `term'_r, replace
		
		// fitted
		qui predict `term'_f, xb equation(`term')
		qui putmata `term'_f = `term'_f, replace
		mata: `term'_f = `term'_f[`start_ind'::`sample_size',.]
		mata: `term'_r = `term'_r[`start_ind'::`sample_size',.]
		
		if (`is_first' == 0) {
		  mata: _res = `term'_r;
		  mata: _fit = `term'_f;
		  local is_first = 1
		}
		else {
		  mata: _res = _res,`term'_r;
		  mata: _fit = _fit,`term'_f;
		}
	  }
	  
	  // Sigma 
	  
	  
	  // Calculate the RMSE
	  mata: RMSE = cal_RMSE(&_res);
	  mata: st_matrix("_RMSE",RMSE);
	  
	  // Calculate R_sq
	  //dis "here"
	  //mata: printf("Fit: %f %f\n",rows(_fit),cols(_fit));
	  //mata: ESS_mat = _fit;// Derive TSS
	  mata: ESS_mat = (_fit)' * (_fit)
	  //mata: printf("Res: %f %f\n",rows(ESS_mat),cols(ESS_mat));
	  mata: mle_ESS = trace(ESS_mat);
	  //mata: printf("mle_ESS: %f \n",mle_ESS);
	  mata: RSS_mat = (_res)' * (_res);
	  mata: st_matrix("res",RSS_mat);
	  mata: mle_RSS = trace(RSS_mat);
	  //mata: printf("mle_RSS: %f \n",mle_RSS);
	  mata: mle_R_sq = mle_ESS / (mle_ESS + mle_RSS);
	  //dis "R^2:"
	  //mata: printf("mle_R_sq: %f \n",mle_R_sq);
	  mata: st_matrix("mle_R_sq_",mle_R_sq);
	  mata: st_matrix("mle_RSS_",mle_RSS);
	  local R_sq = mle_R_sq_[1,1]
    }
	else {
	  local is_success = 0
	  display as error "Error in MLE Estimation"
	  ereturn scalar is_success = 0
	  exit
	}
  restore
  if (`is_success' == 1) {
  local df            = (`sample_size' - `nlag') - 2 * e(k_1)
  
  matrix beta         = e(b)
  //dis "here A"
  matrix Sigma        = e(Sigma)
  //dis "here B"
  matrix sample_Sigma = e(Sigma) / `df'
  //dis "here C"
  local det_Sigma = det(e(Sigma))
  //dis "here D"
  local trace_Sigma = trace(e(Sigma))
  
  ereturn clear
  ereturn matrix beta  = beta
  ereturn matrix Sigma = Sigma
  ereturn matrix res   = res
  ereturn matrix sample_Sigma = sample_Sigma
  
  ereturn scalar df      = `df'
  ereturn scalar Sigma_det = `det_Sigma'
  ereturn scalar RMSE = _RMSE[1,1]
  ereturn scalar Sigma_trace = `trace_Sigma'
  ereturn scalar R_sq = `R_sq'
  ereturn scalar is_success = `is_success'

  }
  /////////////////
  /*0
  ereturn matrix beta    = beta
  ereturn matrix res     = err_mat
  ereturn matrix Sigma = Sigma
  ereturn matrix sample_Sigma = sample_Sigma
  
  ereturn scalar df      = `degree_of_freedom'
  ereturn scalar Sigma_trace = Sigma_trace[1,1]
  ereturn scalar Sigma_det   = Sigma_det[1,1]
  ereturn scalar RMSE = RMSE[1,1]
  ereturn scalar R_sq = R_sq[1,1]
  ereturn scalar is_success = `is_success'
  */
end
  
