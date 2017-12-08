preserve
clear
set obs 30
local num_regressors = 4      // Number of Regressors
local b_times        = 10     // Number of Times to run the simulation (>= 9)
local b_nsample      = _N     // Number of size bootstrap sample (must be >= current sample size)
local is_constant    = 1      // Include a constant term? (Yes = 1, No = 0)
local nlag           = 1      // Nnumber of lags
local d              = 1      // Number of Delays
local initial_obs    = 3      // Index of the initial Observation
local ptrim          = 0.1    // Percentage of trimming
local lv_sig         = 0.05   // Level of Significance

sim_2r_normal `num_regressors' `b_times' `b_nsample' `is_constant' `nlag' `d' `initial_obs' `ptrim' `lv_sig'
restore
