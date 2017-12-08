{smcl}
{* 8apr2017}{...}
{hline}
Help for {hi:TVAR_2r_grid_search (eclass)}{right:(Project 2 - Econ 853 Applied Time Series - Winter 2017)}
{right:Wing Kwong Tsang}
{hline}

{title:Threshold Vector Autoregression with 2 Regimes ( TVAR(2) ): Grid Search}

{title:Syntax}

{p 6 14}{cmd:TVAR_2r_grid_search} {it:varlist} [{cmd:if} {it:exp}] [{cmd:in} {it:range}], 
        {cmd:threshold}({it:varname}) {cmd:ptrim}({it:real}) {cmd:d}({it:integer}) {cmd:nlag}({it:integer}) 
		[{cmd:constant}({it:integer 1}) {cmd:ols}({it:integer 1}) {cmd:criterion}({it:integer 1})
		{cmd:est_table}({it:integer 1})]  

{p}{cmd:TVAR_2r_grid_search} is for use with time-series data.  You must {cmd:tsset} your
data before using {cmd:TVAR_2r_grid_search}; see help {help tsset}.

{p} {it:varlist} and {it:threshold} cannot contain time-series operators.

{title:Description}

{p}{cmd:TVAR_2r_grid_search} implements a grid search of the optimal threshold variable according to a given criterion.
The programme conducts the full estimation of TVAR, including the optimal threshold value and beta coefficients for each regime.
Currently we support 2 criteria: the minimum Residual Sum of Squares and the minimum Determinant of the Sigma matrix.
This programme implements the grid search by estimating the TVAR(2) model with all possible thresholds, and selects
the one according to the criterion. 
The command finds out the optimal threshold and estimates the model with the optimal threshold.

{title:Compulsory Inputs}

{p 0 4}{it: varlist} specifies a list of time-series variables for TVAR(2).

{p 0 4}{bf: threshold} specifies a variable for determining the regime of each observation.

{p 0 4}{bf: ptrim} specifies the percentage of observations to be removed when constructing a 
list of possible threshold values.

{p 0 4}{bf: d} specifies the threshold delay parameter for the model. 
The value must be an integer greater than or equal to 1 and less than or equal to {cmd: nlag}.
For example, a value of 1 means that d the threshold at time t is determined {cmd:d} periods 
ahead of time.

{p 0 4}{bf: nlag} specifies the maximum lag order of the model.
The value must be an integer greater than 1.
For instance, {cmd: nlag}(3) means that the model has lags from order 1 to 3.  

{title:Options}

{p 0 4}{cmd:constant} specifies whether to include a constant in the model.
A value of 0 means excluding a constant; a value of 1 means including a constant.

{p 0 4}{cmd:ols} specifies whether using the OLS or ML estimation. 
A value of 1 means using OLS while 0 means using ML.

{p 0 4}{cmd:criterion} specifies the criterion for the grid search. 
A value of 1 means the minimum root mean square error while 0 means the minimum determinant 
of the Sigma matrix.

{p 0 4}{cmd:est_table} specifies whether to include the estimation table
for TVAR(2) using the optimal threshold. The default value 0 means 
that the programme will not display the table, 1 otherwise.

{title:eReturn}

{bf:Scalars}:

 1. e(num_success)        : number of successful trials
 2. e(minStat)            : the statisitics corresponding to the optimal threshold
 3. e(minStat_thrd)       : the optimal threshold value
 4. e(minStat_thrd_index) : the index of the optimal threshold 

{bf:Matrices}:

Suppose we have n variables, T observations, p lags, c such that 
c is 1 if we have a constant and 0 otherwise, and t possible thresholds.
Dimensions of the matrics are in the parentheses.

 1. e(Stat)               : contains the threshold value, index of the threshold, 
                            the RSS and the determinant of the Sigma matrix of 
							each threshold
                            (t x 5)
 2. e(minStat_r_Sigma)    : the Sigma matrix with the optimal threshold 
                            (n x n)
 3. e(minStat_r1_Sigma)   : the Sigma matrix of the VAR model running with
                            samples from only regime 1
                            (n x n)
 4. e(minStat_r0_Sigma)   : the Sigma matrix of the VAR model running with
                            samples from only regime 2
                            (n x n)
 5. e(minStat_no_r_beta)  : coefficients of the no-Regime model                 
                            (1 x (n*p + c) * n)
 6. e(minStat_r1_beta)    : coefficients of the Regime 1 model                 
                            (1 x (n*p + c) * n)
 7. e(minStat_r0_beta)    : coefficients of the Regime 0 model                 
                            (1 x (n*p + c) * n)
 8. e(minStat_indicator)  : indicator of each observation under the optimal threshold
                            see {help TVAR_2r}.
                            (T x 1)

{title:Associated Files}

Please run the following .ado files before executing the programme {cmd:TVAR_2r_grid_search}

 1. var_ols.ado
 2. var_mle.ado
 3. TVAR_est.ado
 4. grid_search.ado
 5. grid_search_op.ado							
							
{title:Examples}

{bf:Setup}

. {stata "webuse m1gdp": webuse m1gdp}
. {stata "tsset t": tsset t}
. {stata "drop if ln_m1==.": drop if ln_m1==.}
. {stata "gen dln_m1 = d.ln_m1": gen dln_m1 = d.ln_m1}
. {stata "gen dln_gdp = d.ln_gdp": gen dln_gdp = d.ln_gdp}

1. Perform a grid search for the optimal threshold value for TVAR consisting 
   of {it:dln_gdp} and {it:dln_m1} (with constant). 
   {it:dln_m1} is the threshold variable. Trim is set to 15%. 
   Both lag order and delay parameter are set to 1. 
   The estimation is conducted via OLS.
   We also impose restriction on sample period.

. {stata "TVAR_2r_grid_search dln_gdp dln_m1 if t>=tq(1980q1), threshold(dln_m1) ptrim(0.15) d(1) nlag(1) constant(1) ols(1)": TVAR_2r_grid_search dln_gdp dln_m1 if t>=tq(1980q1), threshold(dln_m1) ptrim(0.15) d(1) nlag(1) constant(1) ols(1)}

For other samples, please refer to .do files accompanied with the report.

{title:References}

{p}1. Hansen, B. (1999). Testing for linearity. Journal of Economic Surveys, 13(5):551–576.

{p}2. Hansen, B. E. (1996a). Estimation of TAR Models. Boston College Working Papers in Economics 325.,
      Boston College Department of Economics.

{p}3. Hansen, B. E. (1996b). Inference when a nuisance parameter is not identified under the null hypothesis.
      Econometrica, 64(2):413–430.

{p}4. Lo, M. C. and Zivot, E. (2001). Threshold Cointegration And Nonlinear Adjustment To The Law Of One Price. Macroeconomic Dynamics, 5(04):533–576.

{p}5. Tsay, R. S. (1998). Testing and modeling multivariate threshold models. Journal of the American Statistical Association, 93(443):1188–1202.



