// Final Version
// Can be made quicker by using pointers

//capture mata: mata drop X rt dependent Sigma biased_Sigma biased_Sigma_trace RMSE RSS TSS R_sq // drop mata variables to avoid re definition
//capture mata: mata drop print_matrix() check_X() cons_X() ols_est() cons_err_mat() cons_Sigma() cal_RMSE() cal_TSS()


capture mata: mata drop print_matrix() 
capture mata: mata drop check_X() 
capture mata: mata drop cons_X() 
capture mata: mata drop ols_est() 
capture mata: mata drop cons_err_mat() 
capture mata: mata drop cons_Sigma() 
capture mata: mata drop cal_RMSE() 
capture mata: mata drop cal_TSS()

capture program drop var_ols
program define var_ols, eclass
  version 13
  syntax varlist, nlag(integer) [constant(integer 1) sample_size(integer -1)]
  
  //dis "in OLS"
  local cap_N = _N
  local is_success = 1
  
  //dis "sample size: `sample_size' Vs. `cap_N'"
  if (`sample_size' == -1 | `sample_size' == _N) {
    local sample_size = _N
  }
  else {
    //dis "in here"
    preserve
	qui drop if _n > `sample_size'
  }
  
  local cap_N = _N
  //dis "After sample size: `sample_size' Vs. `cap_N'"
  // Count how many variables in the varlist
  local vars `varlist'
  tokenize `vars'
  local length_varlist = 0
  while "`*'" != "" {
    local current_var `1'
	local length_varlist = `length_varlist' + 1
    macro shift
  }
  
  // Check if df <= 0
  //local degree_of_freedom = `sample_size' - (`length_varlist' * `nlag') - `constant' 
  local degree_of_freedom = `sample_size'
  ereturn scalar df = `degree_of_freedom'
  if (`degree_of_freedom' <= 0) {
    display as error " "
    display as error "Degree of Freedom is 0"
	display as error " "
	local is_success = -1
	ereturn scalar is_success = `is_success' 
	exit
  }
  local num_eqn = `length_varlist'
  qui putmata dependent = (`varlist'), replace
  // Construct matrix X and check validity of X
  mata: X = cons_X(&dependent,`constant', `nlag');
  mata: rt = check_X(&X,rows(X),`length_varlist',`constant');
  mata: st_matrix("_rt",rt);
  local is_success = _rt[1,1]
  
  if (`is_success' != 1) {
    //dis "hello"
    ereturn scalar is_success = `is_success'
    exit
	//dis "after exit"
  }
  //dis "A"
  // Estimation Process
  mata: beta = ols_est(&dependent, &X, `nlag', `constant');
  mata: st_matrix("beta",beta) // Create a stata matrix beta from mata matrix beta 
  
  // construct headers for beta
  local beta_header ""
  foreach dep in `varlist' {
    if (`constant' == 1) {
	  local beta_header `beta_header' "`dep':_cons"
	}
    foreach regressor in `varlist' {
	  forvalues i = 1/`nlag' {
	    local beta_header `beta_header' "`dep':L`i'.`regressor'"
	  }
	}
  }
  matrix colnames beta = `beta_header'
  
  // Construct error matrix
  // col: dependent variable
  mata: err_mat = cons_err_mat(&X, &dependent, &beta,`nlag');
  mata: st_matrix("err_mat",err_mat) // Create a stata matrix beta from mata matrix beta 
  matrix colnames err_mat = `varlist'
  
  // Construct the Sigma matrix and trace
  mata: Sigma = cons_Sigma(&err_mat);
  mata: st_matrix("Sigma",Sigma);
  mata: sample_Sigma = Sigma / (`degree_of_freedom')
  mata: st_matrix("sample_Sigma",sample_Sigma);
  mata: Sigma_trace = trace(Sigma);
  mata: st_matrix("Sigma_trace",Sigma_trace);
  mata: Sigma_det = det(Sigma);
  mata: st_matrix("Sigma_det",Sigma_det);
  matrix colnames Sigma = `varlist'
  matrix rownames Sigma = `varlist'
  
  // Calculate the RMSE
  mata: RMSE = cal_RMSE(&err_mat);
  mata: st_matrix("RMSE",RMSE);

  if (!(`sample_size' == -1 | `sample_size' == _N)) {
    restore
  }
  
  // Calculate the R^2
    // calculate RSS
	mata: RSS = Sigma_trace;
	// Calculate TSS
    mata: TSS = cal_TSS(&dependent,RSS);
  mata: R_sq = 1-(RSS/TSS);
  mata: st_matrix("R_sq",R_sq);
  
  ereturn clear
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
end
// clear mata

mata:

void print_matrix(X) {
  real scalar row_size, col_size;
  row_size = rows((*X));
  col_size = cols((*X));
  for (row = 1; row <= row_size; row++) {
    for (col = 1; col <= col_size; col++) {
	  printf("%3.0f ",(*X)[row,col]);
	}
	printf("\n");
  }
}

real scalar check_X(X, row_size, num_var, is_constant) {
  real scalar rt, num_X_columns, tmp, rank_X;
  rt = 1; // 1 -> no error , 0 -> error
  num_X_columns = num_var + is_constant;
  rank_X = rank(*X);
  //printf("Check: X_col: %f | rank(X): %f | Row(X): %f | Col(X): %f\n",num_X_columns,rank_X, rows(X), cols(X));
  //printf("Rank: %f\n",rank_X);
  correct_rank = cols(*X)
  if (rank_X != correct_rank) {
    rt = 2;
	displayas("err");
	print_matrix(X);
    printf("Design matrix X does not have a full rank (rank(X)=%f but it should be %f)\n",rank_X,correct_rank);
	return (rt)
  }
  return (rt);
}

real matrix cons_X(V, is_constant, nlag) {
  real matrix X, tmp;
  nrows_V = rows(*V);
  ncols_V = cols(*V);
  nrows_X = nrows_V - nlag;
  ncols_X = is_constant + nlag;
  
  col_counter = 0;
  for (V_col = 1; V_col <= ncols_V; V_col++) {
    for (lag = 1; lag <= nlag; lag++) {
	  col_counter++;
	  start_index = nlag - (lag - 1);
	  end_index   = nrows_V - lag;
	  tmp = (*V)[start_index::end_index,V_col];
	  if (col_counter == 1) {
	    X = tmp;
	  }
	  else {
	    X = X,tmp;
	  }
    }
  }
  if (is_constant == 1) {
    tmp = J(rows(X),1,1);
	X = tmp,X;
  }
  //printf("\n\nMatrix X\n");
  //print_matrix(X);
  return (X);
}

real matrix ols_est(V, X, nlag, is_constant) {
  real scalar start_col, row_size, col_size;
  real matrix Y, beta, all_beta;
  V_row_size = rows(*V);
  V_col_size = cols(*V);
  X_row_size = rows(*X);
  X_col_size = cols(*X);
  
  //printf("X: row = %f | col = %f \n",rows(X),cols(X));
  //Y = V[nlag+1..V_row_size,1];
  //printf("Y: row = %f | col = %f \n",rows(Y),cols(Y));
  is_first_beta = 1;
  for (dep = 1; dep <= V_col_size; dep++) {
    Y = (*V)[nlag+1..V_row_size,dep];
	beta = (luinv((*X)' * (*X)) * ((*X)') ) * Y
	//print_matrix(beta);
	if (is_first_beta == 1) {
	  is_first_beta = 0;
	  all_beta = beta';
	}
	else {
	  all_beta = all_beta,beta';
	}
  }  
  return (all_beta);
}

real matrix cons_err_mat(X,V,all_beta,nlag){
  real matrix err_v, err_mat, Y, beta, pred;
  V_num_cols = cols(*V);
  V_num_rows = rows(*V);
  X_num_cols = cols(*X);
  X_num_rows = rows(*X);
  all_beta_start = -1;
  all_beta_end   = -1;
  is_first_err_mat = 0;
  for (dep = 1; dep <= V_num_cols; dep++) {
    Y = (*V)[nlag+1..V_num_rows,dep];
	all_beta_start = (dep - 1) * X_num_cols + 1;
	all_beta_end   = all_beta_start + X_num_cols - 1; 
	beta = ((*all_beta)[1,all_beta_start::all_beta_end])';
	pred = (*X) * (beta);
	err_v = Y - pred;
	
	if (is_first_err_mat == 0) {
	  err_mat = err_v;
	  is_first_err_mat = 1;
	}
	else {
	  err_mat = err_mat,err_v;
	}
  }
  return (err_mat)
}

real matrix cons_Sigma(err_mat) {
  num_cols = cols((*err_mat));
  real matrix Sigma, tmp1, tmp2;
  Sigma = J(num_cols,num_cols,.);
  
  for (a = 1; a <= num_cols; a++) {
    tmp1 = (*err_mat)[.,a];
    for (b = 1; b <= num_cols; b++) {
	  tmp2 = (*err_mat)[.,b];
	  Sigma[a,b] = (tmp1') * tmp2;
	}
  }
  return (Sigma);
}

real scalar cal_RMSE(err_mat) {
  num_cols = cols((*err_mat));
  num_rows = rows((*err_mat));
  num      = num_cols * num_rows;
  sum_err = 0;
  for (row = 1; row <= num_rows; row++) {
    for (col = 1; col <= num_cols; col++) {
	  sum_err = sum_err + (*err_mat)[row,col] * (*err_mat)[row,col];
	}
  }
  return (sqrt(sum_err/num));
}

real scalar cal_TSS(V, RSS) {
  real scalar TSS;
  cols = cols(*V);
  rows = rows(*V);
  all_sum = 0;
  for (i = 1; i <= cols; i++) {
    for (j = 1; j <= rows; j++) {
	  all_sum = all_sum + (*V)[j,i]; 
	}
  }
  _mean = all_sum / (cols*rows);
  all_sum = 0;
  for (i = 1; i <= cols; i++) {
    for (j = 1; j <= rows; j++) {
	  all_sum = all_sum + ((*V)[j,i] - _mean)^2; 
	}
  }
  return (all_sum);
}
end
