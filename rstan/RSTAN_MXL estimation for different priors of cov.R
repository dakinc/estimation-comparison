#####-----#####-----#####-----#####
## D.Akinc - June 2017
## R code to execute STAN estimation: Hierarchical Bayes Estimation
## 
## The "fullData.csv" follows the long-format style used in "R-package:mlogit"
## see Croissant, Y. (2012). Estimation of multinomial logit models in R: 
## The mlogit Packages. R package version 0.2-2. 
## URL: http:// cran.r-project.org/web/packages/mlogit/vignettes/mlogit.pdf.
## e.g.
## [id | alt | Att11 | Att12 | Att21 | Att22 | Att31 | Att32 | choice {0,1}]
##
## Remark that the "fullData.csv" contains all data for 10 repetitions in long-format
## index1 = the row no. of the 1st observation in each repetition
## index2 = the row no. of the last observation in each repetition
## index3 = id of the repetition

#####-----#####-----#####-----#####
mxlstan = function(index1,index2,index3){
Sys.time()->start;
setwd("/rstan")

## Required libraries:
library(data.table)
library(mlogit)
library(rstan)
library(doParallel)
library(parallel)
  
#numCores<-(detectCores()/2)

#Dataset setup
nalts=3 #num of alternatives
nsets=18 #num of choice sets
nres=200 #num of respodents
df=6 #num of model parameters
nrep=10 #num of repetition

#MCMC Settings
nsample = 15000
nburnin = 35000
nitit = nsample + nburnin

#Storage for output
RMSEmatrix=NULL
estbetafinal=NULL
estmufinal=NULL
estsigmafinal=NULL
  
#Importing data
fullData=read.table("/fullDATA.csv",header=TRUE, sep=";")

#####-----#####-----#####-----#####
## Read in the choice data and prepare for estimation
#####-----#####-----#####-----#####
kkk=index1:index2
chdat <- fullData[kkk,]
dat <- mlogit.data(data = chdat, id = "id", choice = "choice", 
				   shape = "long", alt.var = "alt", varying=4:9)
dat <- data.table(dat)
dat[, alt := as.integer(alt)]
Formula_m = formula(~ Att11 + Att12 + Att21 + Att22 + Att31 + Att32 - 1)
  
dat_choice = dat[choice == 1L]
S = nrow(dat_choice)
P = length(all.vars(update(formula_m, NULL ~ .)))
K = 3L
N = uniqueN(dat_choice[['id']])
D = 1L
dem = matrix(1, N)
nu = 2L
  
dat_u <- as.data.frame(model.matrix(~ Att11 + Att12 + Att21 
									+ Att22 + Att31 + Att32 - 1, data = dat))
x_arr = array(data = unname(model.matrix(formula_m, data = dat_u)),dim = c(K, S, P))
x_arr = aperm(x_arr, c(2, 1, 3))
  
#####-----#####-----#####-----#####
## Input data for RSTAN 
#####-----#####-----#####-----#####
dat_rstn =
    list(S = S, P = P, D = D, K = K, N = N, # dimensions
         y = dat_choice[['alt']],   # choices = dependent variable
         id = dat_choice[['id']],   # repondents ids
         x = x_arr,                 # design matrix
         d = dem)                   # group level predictors - not used, therefore 1
		 #nu = nu) 					# degrees of freedom for the half-t prior distribution
#####-----#####-----#####-----#####  
## RSTAN model:
## (1) Call 'MXLestIW.stan' file for inverse Wishart prior
## (2) Call 'MXLestSS.stan' file for Separation  Strategy with the setting of 
##	   Barnard, et al.(2000)
## (3) Call 'MXLestSS_STAN.stan' file for Separation  Strategy with the setting of 
##	   STAN manual(2016)
## (4) Call 'MXLestSIW.stan' file for scaled inverse Wishart prior
## (5) Call 'MXLestHALFT.stan' file for Huang Half-t prior
#####-----#####-----#####-----#####
res_stan = stan(file = 'MXLestIW.stan',
                 data = dat_rstn,
                 chains = 2L, cores = 2L, iter = nitit, warmup = nburnin,
                 verbose = TRUE, pars = c('mu', 'Sigma', 'beta', 'lp__'))
#print(res_stan, pars = "mu", digits_summary = 3)
#stan_trace(res_stan, pars = "mu", inc_warmup = TRUE)
#print(res_stan,par = c("mu", "Sigma", "lp__"), digits_summary = 4)
  
#####-----#####-----#####-----#####
## Results from Estimation
#####-----#####-----#####-----#####
eststan <- extract(res_stan, pars = c("mu", "Sigma","beta"), permuted=TRUE)

#True Population parameters used in simulating choice data  
truebeta  <- read.table("/betatrue.csv",header=TRUE, sep=";")
truebeta <- t(truebeta)
truemu  <-0.5*c(-1,0,-1,0,-1,0)
sigmatrue=0.25*diag(df)
nlong = (df*(df+1))/2;
sigmatruelong = matrix(0,nlong,1)
j2 = 0;
for (s in 1:df){
	j1 = j2 + 1
    j2 = j1 + (df - s)
    sigmatruelong[j1:j2,] = sigmatrue[s,s:df]}  
#####-----#####-----#####-----#####
## Computation of RMSEs
#####-----#####-----#####-----#####
#RMSE_BETA:
estbeta = colMeans(eststan$beta)
estbeta=t(estbeta)
betasum = 0
for (m in 1:nres){
    betasum = betasum + (t(estbeta[,m]-truebeta[,m])
			%*%(estbeta[,m]-truebeta[,m]))}
RMSEbeta=sqrt(betasum/(nres*df))
#RMSE_MU
estmu=colMeans(eststan$mu)
RMSEmu=sqrt((t(estmu - truemu)%*%(estmu - truemu))/df)
#RMSE_SIGMA
estsigma = colMeans(eststan$Sigma)
sigmameanlong = matrix(0,nlong,1)
j2 = 0;
for (t in 1:df){
    j1 = j2 + 1
    j2 = j1 + (df - t)
    sigmameanlong[j1:j2,] = estsigma[t:df,t]}
RMSEsigma=sqrt((t(sigmameanlong - sigmatruelong)%*%(sigmameanlong - sigmatruelong))/nlong)
RMSEsim <- c(RMSEbeta, RMSEmu, RMSEsigma)
  
#Final Estimates Matrices
estbetafinal=rbind(estbetafinal,estbeta)
estmufinal=rbind(estmufinal,estmu)
estsigmafinal=rbind(estsigmafinal,sigmameanlong)
RMSEmatrix=rbind(RMSEmatrix,RMSEsim)
print(Sys.time()-start);
  
write.table(estbetafinal, file = paste("/estbeta_","rep",index3,".csv",sep=""), row.names = FALSE,sep=";")
write.table(estmufinal, file = paste("/estmu_","rep",index3,".csv",sep=""), row.names = FALSE,sep=";")
write.table(estsigmafinal, file = paste("/estsigma_","rep",index3,".csv",sep=""), row.names = FALSE,sep=";")
write.table(RMSEmatrix, file = paste("/RMSE_","rep",index3,".csv",sep=""), row.names = FALSE,sep=";")
cat("Estimation Procedure is finalized")
}
mxlstan(1,10800,01)
mxlstan(10801,21600,02)
mxlstan(21601,32400,03)
mxlstan(32401,43200,04)
mxlstan(43201,54000,05)
mxlstan(54001,64800,06)
mxlstan(64801,75600,07)
mxlstan(75601,86400,08)
mxlstan(86401,97200,09)
mxlstan(97201,108000,10)
