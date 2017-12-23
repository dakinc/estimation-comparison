#####-----#####-----#####-----#####
## D.Akinc - June 2017
## R code for the Hierarchical Bayes Estimation of the Mixed Logit Model
## by using the Inverse Wisharts prior for the covariance matrix
## Here the function "rhierMnlRwMixture" of bayesm package is used 
## for estimating the mixed logit model.
## The computation time of the function rhierMnlRwMixture is significantly faster
## than the "MXL - MCMC estimation - IW", it is preferable to use the bayesm package 
## when using IW priors for the covariance matrix.
## 
## The file "xmat.csv" follows the long-format style
## e.g.
## [Att11 | Att12 | Att21 | Att22 | Att31 | Att32 ]
## The file "ymat.csv" contains the choices of each respondent for each choice set.
##

library(MASS)
# install.packages("bayesm")
require(bayesm)

rm(list=ls())
Sys.time()->start;

# Dataset setup
nalts = 3 #num of alternatives
nsets = 18 #num of choice sets
nres = 200 #num of respodents
df = 6 #num of model parameters
pkg = 1

# MCMC Settings
nsample = 30000
nburnin = 70000
R = nsample + nburnin
keep = 1
Mcmc1=list(R=R,keep=keep)
ind= R - nsample + 1
totalit=R/keep

# Hyperparameters for Inverse Wishart
# The hyperparameters of IW priors used in the Simulation Study, respectively:
# (1) IW1 --> nu = df & V = diag(df)
# (2) IW2 --> nu = df & V = nu * diag(df)
# (3) IW3 --> nu = (df+1) & V = nu * diag(df)
# (4) IW4 --> nu = (df+3) & V = nu * diag(df)
# (5) IW5 --> nu = (df+4) & V = nu * diag(df)
# (6) IW6 --> nu = 0.5 * (df * (df+1)) & V = 0.1 * nu *diag(df)
nu = df  #nu = p
V = nu * diag(df) #T = v*I
Prior1=list(ncomp=1,nu=nu,V=V)

# Importing data
xmat <- read.table("/xmat.csv",header=TRUE,sep=";")
ymat <- read.table("/ymat.csv",header=TRUE,sep=";")
xmat=t(t(xmat))
ymat=t(t(ymat))
i=1 
xind=((i-1)*nsets*nalts*nres+1):(i*nsets*nalts*nres) 
yind=((i-1)*nres+1):(i*nres)
x=xmat[xind,]
y=ymat[yind,]
yvec=as.vector(t(y))
yvec=t(t(yvec))

# R code used for converting data into the data.frame format reqired for bayesm package
source("/userin.R",local=TRUE)
testData=userin(pkg,x,yvec,nalts,nsets,nres,df)

# True Population parameters used in simulating choice data  
truebeta <- read.table("/betatrue.csv",header=T, sep=";")
truebeta = t(truebeta)
truemu  = 2.0*c(-1.0,0.0,-1.0,0.0,-1.0,0.0)
sigmatrue = 1.0*diag(df)
nlong = (df*(df+1))/2;
sigmatruelong = matrix(0,nlong,1)
j2 = 0;
for (s in 1:df){
  j1 = j2 + 1
  j2 = j1 + (df - s)
  sigmatruelong[j1:j2,] = sigmatrue[s:df,s]}

#####-----#####-----#####-----#####
## Estimation 
#####-----#####-----#####-----#####
estout = rhierMnlRwMixture(Data=testData,Prior=Prior1,Mcmc=Mcmc1)

#####-----#####-----#####-----#####
## Results from Estimation
#####-----#####-----#####-----#####
estbeta <- apply(estout$betadraw[,,ind:totalit], c(1,2), mean)
estbeta = t(estbeta)

mudraw=estout$nmix$compdraw
postmu=NULL
for (k in ind:totalit){
  postmu=rbind(postmu,mudraw[[k]][[1]]$mu)
}
estmu = colMeans(postmu)

postsigma=matrix(0,df,df)
for (l in ind:totalit){
  postsigma=postsigma+crossprod(backsolve(mudraw[[l]][[1]]$rooti,diag(df)))
}
estsigma=postsigma/nsample

#####-----#####-----#####-----#####
## Computation of RMSEs
#####-----#####-----#####-----#####
# RMSE_BETA:
betasum = 0
for (m in 1:nres){
	betasum = betasum + (t(estbeta[,m]-truebeta[,m])
			%*%(estbeta[,m]-truebeta[,m]))}
RMSEbeta = sqrt(betasum/(nres*df))

# RMSE_MU
RMSEmu = sqrt((t(estmu - truemu)%*%(estmu - truemu))/df)

# RMSE_SIGMA
sigmameanlong = matrix(0,nlong,1)
j2 = 0;
for (t in 1:df){
	j1 = j2 + 1
    j2 = j1 + (df - t)
    sigmameanlong[j1:j2,] = estsigma[t:df,t]}
RMSEsigma = sqrt((t(sigmameanlong - sigmatruelong)%*%(sigmameanlong - sigmatruelong))/nlong)
RMSEmatrix = cbind(RMSEbeta,RMSEmu,RMSEsigma)

print(Sys.time()-start);

# Exportin the Estimates and the RMSE values
# write.table(estbeta, file = "/estbeta.csv", row.names = FALSE,sep=";")
# write.table(estmu, file = "/estmu.csv", row.names = FALSE,sep=";")
# write.table(estsigma, file = "/estsigma.csv", row.names = FALSE,sep=";")
# write.table(RMSEmatrix, file = "/RMSE.csv", row.names = FALSE,sep=";")
cat("Estimation Procedure is finalized")