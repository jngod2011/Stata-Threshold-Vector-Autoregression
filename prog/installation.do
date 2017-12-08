// This installation file is for Windows users.
// Please also ensure that the current working directory contains all files
// and the file/directory structure is preserved

// Package Install do file

// install estimation files
  do "estimation\TVAR_2r.ado"
  do "estimation\TVAR_est.ado"
  do "estimation\var_mle.ado"
  do "estimation\var_ols.ado"


// install grid search files
  do "grid_search\grid_search_op.ado"
  do "grid_search\grid_search.ado"
  
// install Monte Carlo Simiulation files
  do "monte\bootstrap.ado"
  do "monte\normal_simulation.ado"
