In R folder,

Simulation data generator.R
 - Generating simulation data based on our simulation setting
 - Please refer censoring c.0 ~ c.7 at bottom of file.
 - For example, c.5 is censor parameter that makes 50% censoring proportion.
 
Weight generator.R
  - Generating function of IPCW weight (Default is using survival function of WKM function of emplik package)
  
is_est.R
  - Our suggested method to estimate beta using nleqslv function
  - Paper (Doctor ver.) equation (4)
  
  
is_optim_est.R	Revised on 20-06-08	3 days ago
  - Our suggested method to estimate beta using optim function
  - Paper (Doctor ver.) equation (5)

rq_est.R
  - suggested method by Kim et al.(2012) to estimate beta using nleqslv function 
  - Censored quantile regression for residual lifetimes (Kim et al. 2012) equation (9)

rq_optim_est.
  - suggested method by Kim et al.(2012) to estimate beta using optim function
  - Censored quantile regression for residual lifetimes (Kim et al. 2012) equation (8)

ndata2.R
 - simulation example
 - setting : median quantile estimation from dataset with various censoring proportion (0% ~ 70%)
 
Reference papers.
Regression on Quantile Residual life (Jung et al. 2009).pdf
  - paper about simulation settings

Censored quantile regression for residual lifetimes (Kim et al. 2012).pdf
  - paper about unsmoothed estimating equation for quantile regression
  
Paper (Doctor ver.).pdf
  - Our paper about smoothed estimating equation for quantile regression
  
 
