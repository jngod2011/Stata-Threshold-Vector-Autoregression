{smcl}
{* 8apr2017}{...}
{hline}
Help for {hi:TVAR_2r (eclass)}{right:(Econ 853 Applied Time Series - Winter 2017)}
{right:Wing Kwong Tsang}
{hline}

{title:Threshold Vector Autoregression with 2 Regimes ( TVAR(2) ): Estimation}

{p 6 14}{cmd:TVAR_2r} {it:varlist} [{cmd:if} {it:exp}] [{cmd:in} {it:range}], {cmd:indicator}({it:varname}) {cmd:nlag}({it:integer}) {cmd:d}({it:integer}) [{cmd:constant}({it:integer 1}) {cmd:ols}({it:integer 1}) {cmd:no_regime}({it:integer 1})]

{p}{cmd:TVAR_2r} is for use with time-series data.  You must {cmd:tsset} your
data before using {cmd:TVAR_2r}; see help {help tsset}.

{p} {it:varlist} and {it:indicator} cannot contain time-series operators.

{title:Description}

{cmd:TVAR_2r} estimates a TVAR(2) (with the first regime indexed as 0 and the second regime indexed as 1)
using the Ordinary Least Square (OLS) Estimation or the Maximum Likelihood (ML) Estimation. 
Given 1 or more time series, the indicator variable, 
the number of lags and the delay parameter, one can obtain the root mean square error (RMSE), 
R^2 and coefficients for the Regime 0 model, Regime 1 model and the no-Regime Model.

{p}One can also specify whether the model has a constant term, whether using OLS or ML for estimation 
and whether to display the result of the no-Regime model.

The user {title:must provide an indicator variable} to indicate which Regime the observation belongs to. 
This allows users to specify the model using their own criteria.

{title:Compulsory Inputs}

{p 0 4}{it: varlist} specifies a list of time-series variables for TVAR(2).

{p 0 4}{cmd: indicator} specifies a variable that signifies the regime of each observation
(i.e. a threshold variable).
A value of 0 signifies that the observation belongs to Regime 0 while a value of 1 means 
that the observation belongs to Regime 1.

{p 0 4}{cmd: nlag} specifies the maximum lag order of the model.
The value must be an integer greater than 1.
For instance, {cmd: nlag}(3) means that the model has lags from order 1 to 3.  

{p 0 4}{cmd: d} specifies the threshold delay parameter for the model. 
The value must be an integer greater than or equal to 1 and less than or equal to {cmd: nlag}.
For example, a value of 1 means that d the threshold at time t is determined {cmd:d} periods 
ahead of time.

{title:Options}

{p 0 4}{cmd:constant} specifies whether to include a constant in the model.
A value of 0 means excluding a constant; a value of 1 means including a constant.

{p 0 4}{cmd:ols} specifies whether using the OLS or ML estimation. 
A value of 1 means using OLS while 0 means using ML.

{p 0 4}{cmd:no_regime} specifies whether to show the estimation result of the 
no-Regime model. A value of 1 means showing the estimation results of the 
no-Regime model; a value of 0 means not showing.

{title:eReturn}

{bf:Scalars}:

 1. e(r0_is_success)  : 1 if estimation of Regime 0 is successful; 0 otherwise.
 2. e(r0_df)          : degree of freedom of Regime 0 estimation.
 3. e(r0_size)        : sample size of Regime 0
 4. e(r0_trace)       : trace of the Sigma matrix of Regime 0 Estimation.
 5. e(r0_det)         : determinant of the Sigma matrix of Regime 0 Estimation.
 6. e(r0_RMSE)        : the root mean square error (RMSE) of Regime 0 estimation.
 7. e(r0_R_sq)        : the R^2 of the Regime 0 model.
 
 8. e(r1_is_success)  : 1 if estimation of Regime 1 is successful; 0 otherwise.
 9. e(r1_df)          : degree of freedom of Regime 1 estimation.
10. e(r1_size)        : sample size of Regime 1
11. e(r1_trace)       : trace of the Sigma matrix of Regime 1 Estimation.
12. e(r1_det)         : determinant of the Sigma matrix of Regime 1 Estimation.
13. e(r1_RMSE)        : the root mean square error (RMSE) of Regime 1 estimation.
14. e(r1_R_sq)        : the R^2 of the Regime 1 model.

15. e(r_Sigma_det)    : the determinant of the overall Sigma matrix 
16. e(r_Sigma_trace)  : the trace of the overall Sigma matrix

17. e(no_r_is_success): 1 if estimation of Regime 1 is successful; 0 otherwise.
18. e(no_r_df)        : degree of freedom of Regime 1 estimation.
19. e(no_r_trace)     : trace of the Sigma matrix of Regime 1 Estimation.
20. e(no_r_det)       : determinant of the Sigma matrix of Regime 1 Estimation.
21. e(no_r_R_sq)      : the root mean square error (RMSE) of Regime 1 estimation.
22. e(no_r_RMSE)      : the R^2 of the Regime 1 model.

{bf:Matrices}:

Suppose we have n variables, T observations, p lags, and c such that 
c is 1 if we have a constant and 0 otherwise.
Dimensions of the matrics are in the parentheses.

 1. e(no_r_beta)      :  coefficients of the no-Regime model                 
                         (1 x (n*p + c) * n)
 2. e(no_r_res)       :  residuals of the no-Regime model                    
                         (T x n)
 3. e(no_r_Sigma)     :  the Sigma matrix of the no-regime model             
                         (n x n)
 4. e(r_Sigma)        :  the overall Sigma matrix of the models with regimes 
                         (n x n)
 5. e(r1_beta)        :  coefficients of the Regime 1 model                  
                         (1 x (n*p + c) * n)
 6. e(r1_res)         :  residuals of the Regime 1 model                     
                         (T x n)
 7. e(r1_Sigma)       :  the Sigma matrix of the Regime 1 model              
                         (n x n)
 8. e(r0_beta)        :  coefficients of the Regime 0 model                  
                         (1 x (n*p + c) * n)
 9. e(r0_res)         :  residuals of the Regime 0 model                     
                         (T x n)
10. e(r0_Sigma)       :  the Sigma matrix of the Regime 0 model              
                         (n x n)

{title:Examples}

Please install associated files before running examples below.

{it:{bf:Setup}}

. {stata "webuse set http://www.stata.com/users/dschenck/":webuse set http://www.stata.com/users/dschenck/}
. {stata "webuse usmacro.dta"     :webuse usmacro.dta}

Assume that observations with an odd index belongs to Regime 1. 
Other observations belong to Regime 0.

. {stata "gen ind = 0"     :gen ind = 0}
. {stata "replace ind = 1 if mod(_n,2) == 1"     :replace ind = 1 if mod(_n,2) == 1}

Fit TVAR(2) by OLS with constant, 1 lag (the default), 0 delay (the default) and dln_m1 being the threshold variable 
from the fifth observation to the 200th observation.

. {stata "TVAR_2r unrate inflation in 5/200, indicator(ind) nlag(1) d(0) constant(1) ols(1)":TVAR_2r unrate inflation in 5/200, indicator(ind) nlag(1) d(0) constant(1) ols(1)}

{title:Associated Files}

Please run the following .ado files before executing the programme {cmd:TVAR_2r_grid_search}

 1. var_ols.ado
 2. var_mle.ado
 3. TVAR_est.ado

{title:References}

{p}1. Hansen, B. (1999). Testing for linearity. Journal of Economic Surveys, 13(5):551–576.

{p}2. Hansen, B. E. (1996a). Estimation of TAR Models. Boston College Working Papers in Economics 325.,
      Boston College Department of Economics.

{p}3. Hansen, B. E. (1996b). Inference when a nuisance parameter is not identified under the null hypothesis.
      Econometrica, 64(2):413–430.

{p}4. Lo, M. C. and Zivot, E. (2001). Threshold Cointegration And Nonlinear Adjustment To The Law Of One Price. Macroeconomic Dynamics, 5(04):533–576.

{p}5. Tsay, R. S. (1998). Testing and modeling multivariate threshold models. Journal of the American Statistical Association, 93(443):1188–1202.

