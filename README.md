## The Threshold Vector-autogression Program Package 

This respository contains a Stata Threshold Vector Auto-regression program package.

## List of Programs:
The .zip file contains:

|     | File Name                                  | Description                                                                                                      |
|-----|--------------------------------------------|------------------------------------------------------------------------------------------------------------------|
| 1.  | prog/installation.do                       | Initialisation of the environment for the program package                                                        |
|     |                                            |                                                                                                                  |
| 2.  | prog/estimation/TVAR_2r.ado                | The main program to estimate a 2-regime TVAR                                                                     |
| 3.  | prog/estimation/TVAR_2r.sthlp              | The help file for TVAR_2r                                                                                        |
| 4.  | prog/estimation/TVAR_est.ado               | A helper program for estimating a TVAR model                                                                     |
| 5.  | prog/estimation/var_mle.ado                | A helper program to estimate the TVAR using the maximum likelihood                                               |
| 6.  | prog/estimation/var_ols.ado                | A helper program to estimate the TVAR using the ordinary least square                                            |
|     |                                            |                                                                                                                  |
| 7.  | prog/grid_search/grid_search.ado           | The main (wrapper) program to implement grid searches for a TVAR model                                           |
| 8.  | prog/grid_search/grid_search_op.ado        | A helper program for a TVAR model's grid searches                                                                |
| 9.  | prog/grid_search/TVAR_2r_grid_search.sthlp | The help file for TVAR_2r_grid_search                                                                            |
|     |                                            |                                                                                                                  |
| 10. | prog/monte/bootstrap.ado                   | A program to implement bootstrapping re-sampling process                                                         |
| 11. | prog/monte/normal_simulation.ado           | A program to implement a Monte Carlo Simulation for TVAR using random variables of uniform[0,1] and Normal(0,1)  |
| 12. | prog/monte/test_normal_simulation.ado      | A wrapper program for implementing the above Monte Carlo simulation                                              |

## Examples:

Please see the help files in prog/help for how the package works.
