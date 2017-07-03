#####-----#####-----#####-----#####
## D.Akinc - June 2017
## R code for the Hierarchical Bayes Estimation of the Mixed Logit Model
## by using the Scaled Inverse Wisharts prior for the covariance matrix
## 
## The file "xmat.csv" follows the long-format style
## e.g.
## [Att11 | Att12 | Att21 | Att22 | Att31 | Att32 ]
## The file "ymat.csv" contains the choices of each respondent for each choice set.
##
## Remark that both "xmat.csv" and "ymat.csv" contain all data 
## for 10 repetitions in long-format.

#Required libraries for parallel computation:
library(foreach)
library(doParallel)
library(parallel)

#Parameters for parallel computation
numCores<-(detectCores()/2)
cl<-makeCluster(numCores)
registerDoParallel(cl)
Sys.time()->start;

#Dataset setup
nalts=3 #num of alternatives
nsets=18 #num of choice sets
nres=200 #num of respodents
df=6 #num of model parameters
nrep=10 #num of repetition
i<-1:nrep

#MCMC Settings
nsample = 30000
nburnin = 70000
nskip = 1

#Hyperparameters for Inverse Wishart
nu = df  #v0 = p
V = nu * diag(df)   

#Importing data
xmat=read.table("/xmat.csv",header=TRUE,sep=";")
ymat=read.table("/ymat.csv",header=TRUE,sep=";")
xmat=t(t(xmat))
ymat=t(t(ymat))

#####-----#####-----#####-----#####
## Estimation 
#####-----#####-----#####-----#####
comb <- function(x, ...) {
  lapply(seq_along(x),function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))}

AllOutputList<-foreach(i, .combine='comb',.packages = c("MCMCpack", "MASS"),
				.multicombine=TRUE,.init=list(list(), list(), list(), list()))%dopar%{
	tryCatch({
#####-----#####-----#####-----#####
## Functions for MCMC estimation
#####-----#####-----#####-----#####
## Function for computing the loglikelihood of MNL model
loglike = function(nalts,nsets,nres,betan,y,X){
logl = matrix(0,1,nres)
for (r in 1:nres){
	tt = ((r-1)*nalts*nsets+1):(r*nalts*nsets)
	Xsub = X[tt,]
	row = nrow(Xsub)
	nset = row/nalts
	logl[,r] = 0
		for (s in 1:nset){
		xx = Xsub[((s-1)*nalts+1):(s*nalts),]
		#m = max(xx[1,]%*%betan[,r], xx[2,]%*%betan[,r])
		m = max(xx[1,]%*%betan[,r], xx[2,]%*%betan[,r], xx[3,]%*%betan[,r])
		u = (xx%*%betan[,r])-m
		p = exp(u) / sum(exp(u))
		logl[,r] = logl[,r] + log(p[y[s,r],])}}
return(logl)}

## Function for the conditional posterior of MU_BETA
nextmu = function(nres,betan,Sigma){
lb = t(chol(Sigma))
ndf = ncol(Sigma)
zr = matrix(rnorm(ndf),nrow=ndf,ncol=1)
munew = rowMeans(betan)+((lb%*%zr)/sqrt(nres))
return(munew)}

## Function for the conditional posterior of BETA_N
nextbeta = function(nres,betan,rho,mu,Sigma,y,X,logLold){
#Metropolis-Hasting Step:
lb = t(chol(Sigma))
ndf = ncol(Sigma)
zr = matrix(rnorm(df*nres),nrow=df, ncol=nres)
betanew = betan + sqrt(rho)*(lb%*%zr)
bn = betanew - matrix(mu,nrow=ndf,ncol=nres)
bo = betan - matrix(mu,nrow=ndf,ncol=nres)
logLnew = loglike(nalts,nsets,nres,betanew,y,X)
varmat1 = bn * (solve(Sigma)%*%bn)
varmat2 = bo * (solve(Sigma)%*%bo)
r = exp(logLnew - logLold - 0.5 * (colSums(varmat1) - colSums(varmat2)))
ind = (matrix(runif(nres),nrow=1,ncol=nres)<=r)
meanind = rowSums(ind)/nres
newrho = rho - 0.001*(meanind<0.3)+0.001*(meanind>0.3)
logLnew = logLnew*ind+logLold*(1-ind)
bnew = (betanew*1*matrix(rep(ind,ndf),nrow=ndf,ncol=nres,byrow=T))
	 +(betan*1*matrix(rep((1-ind),ndf),nrow=ndf,ncol=nres,byrow=T))
betaoutlist = list(bnew=bnew,logLnew=logLnew,newrho=newrho,meanind=meanind)
return(betaoutlist)}

#------------------------------------------------------------------------#
#Based on the chosen prior for covariance matrix, the function for the 
#conditional posterior  of SIGMA_BETA (named "nextSigma") and the Gibbs
#Sampling function for updating the all parameters (named "rgibbs") change.
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
###(II.) Under the Scaled Inverse Wishart prior
## Function for the conditional posterior of SIGMA_BETA
nextSigma = function(nres,df,nu,V,betan,mu,delta){
S = (betan - matrix(mu,nrow=df,ncol=nres))
Vnew = diag(1/delta)*(S%*%t(S))*diag(1/delta) + V
Phi = riwish(nu+nres, Vnew)
Sigmanew = diag(delta)*Phi*diag(delta)
return(Sigmanew) }

## Function for the conditional posterior of delta_p
nextdelta = function(df,delta,scaledelta,logdeltaold){
#Metropolis-Hasting Step:
repeat {
	tdelta = (rt(df,df=3)*scaledelta) + delta
    if(sum(tdelta<=1)==0) { break } }
deltanew = log(tdelta)
logdeltanew = log(dt((deltanew-delta)/scaledelta,df=3)*(1/scaledelta))
r = exp(logdeltanew - logdeltaold)
ind = (matrix(runif(df),nrow=1,ncol=df)<=r)
inddelta = rowSums(ind)/df
logdeltanew = logdeltanew*ind+logdeltaold*(1-ind)
dnew = (deltanew*1*matrix(rep(ind,df),nrow=1,ncol=df))+(delta*1*matrix(rep((1-ind),df),nrow=1,ncol=df))
deltaoutlist = list(dnew=dnew,logdeltanew=logdeltanew,inddelta=inddelta)
return(deltaoutlist) }

## Function for updating all parameters by using the Gibbs Sampling
rgibbs = function (nalts,nsets,df,nres,nsample,nburnin,nu,V,y,X){
#Storage for all draws 
mudraw = matrix(0,df,nsample)
Sigmadraw = matrix(0,(df*nsample),df)
betamean = matrix(0,df,nres)
#Starting values:
mu = matrix(((runif(1)*4)-2),df,1)
betao = matrix(mu,nrow=df,ncol=nres)
Sigma = (runif(1)*6)*diag(df)
dsig = 0.1
deltao = matrix(exp(rnorm(df,0,dsig)),nrow=1, ncol=df)
mudeltao = matrix(0,nrow=1,ncol=df)
rho = 0.01
scaledelta = 1.0
suminddelta = 0
logLold = loglike(nalts,nsets,nres,betao,y,X)
logdeltaold = log(dt((deltao-mudeltao)/scaledelta,df=3)*(1/scaledelta))
ndraw = nburnin + nsample*nskip
for (r in 1:ndraw){
    newbetout = nextbeta(nres,betao,rho,mu,Sigma,y,X,logLold)
    betao = newbetout$bnew
    logLold = newbetout$logLnew
    rho = newbetout$newrho
    mu = nextmu(nres,betao,Sigma) 
    Sigma = nextSigma(nres,df,nu,V,betao,mu,deltao)
    newdeltaout = nextdelta(df,deltao,scaledelta,logdeltaold)
    deltao = newdeltaout$dnew
    logdeltaold = newdeltaout$logdeltanew
    suminddelta = suminddelta + newdeltaout$inddelta
    nk = (r-nburnin)/nskip
    if(r>nskip & floor(nk)==nk & nk>0){
      mudraw[,nk]=mu
      nnk = ((nk-1)*df+1):(nk*df)
      Sigmadraw[nnk,]=Sigma
      betamean = betamean+betao } }
meaninddelta = suminddelta / ndraw
betamean = betamean / nsample
mumean = rowMeans(mudraw)
sigmamean = matrix(0,df,df)
for(i in 1:nsample){
	k=((i-1)*df+1):(i*df)
    sigmamean = sigmamean + Sigmadraw[k,]}
sigmamean = sigmamean / nsample
return (list(betamean=betamean,mumean=mumean,sigmamean=sigmamean))}
#------------------------------------------------------------------------#
#####-----#####-----#####-----#####
## Results from Estimation
#####-----#####-----#####-----#####
kkk=((i-1)*nres+1):(i*nres)
lll=((i-1)*nsets*nalts*nres+1):(i*nsets*nalts*nres)
x=xmat[lll,]
y=ymat[kkk,]; y=t(y)

#True Population parameters used in simulating choice data  
truebeta <- read.table("betatrue.csv",header=T, sep=";")
truebeta = t(truebeta)
truemu  <-c(-0.5,0.0,-0.5,0.0,-0.5,0.0)
sigmatrue= 0.25*diag(df)
nlong = (df*(df+1))/2;
sigmatruelong = matrix(0,nlong,1)
j2 = 0;
for (s in 1:df){
	j1 = j2 + 1
	j2 = j1 + (df - s)
	sigmatruelong[j1:j2,] = sigmatrue[s:df,s]}

estout=rgibbs(nalts,nsets,df,nres,nsample,nburnin,nu,V,y,x)

#####-----#####-----#####-----#####
## Computation of RMSEs
#####-----#####-----#####-----#####
#RMSE_BETA:
estbeta = estout$betamean
betasum = 0
for (m in 1:nres){
	betasum = betasum + (t(estbeta[,m]-truebeta[,m])
			%*%(estbeta[,m]-truebeta[,m]))}
RMSEbeta=sqrt(betasum/(nres*df))
#RMSE_MU
estmu=estout$mumean
RMSEmu=sqrt((t(estmu - truemu)%*%(estmu - truemu))/df)
#RMSE_SIGMA
estsigma = estout$sigmamean
sigmameanlong = matrix(0,nlong,1)
j2 = 0;
for (t in 1:df){
	j1 = j2 + 1
    j2 = j1 + (df - t)
    sigmameanlong[j1:j2,] = estsigma[t:df,t]}
RMSEsigma=sqrt((t(sigmameanlong - sigmatruelong)%*%(sigmameanlong - sigmatruelong))/nlong)
RMSEsim=cbind(RMSEbeta,RMSEmu,RMSEsigma)

#Final Estimates Matrices
list(RMSEsim,estbeta,estmu,sigmameanlong)
	}, error=function(e){cat("ERROR :",i,conditionMessage(e), "\n")})}
print(Sys.time()-start);
stopCluster(cl)
RMSEmatrix=do.call(rbind,AllOutputList[[1]])
estbetafinal=do.call(rbind,AllOutputList[[2]])
estmufinal=do.call(rbind,AllOutputList[[3]])
estsigmafinal=do.call(rbind,AllOutputList[[4]])

write.table(estbetafinal, file = "/estbeta.csv", row.names = FALSE,sep=";")
write.table(estmufinal, file = "/estmu.csv", row.names = FALSE,sep=";")
write.table(estsigmafinal, file = "/estsigma.csv", row.names = FALSE,sep=";")
write.table(RMSEmatrix, file = "/RMSE.csv", row.names = FALSE,sep=";")
cat("Estimation Procedure is finalized")

