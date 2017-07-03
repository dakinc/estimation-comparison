#####-----#####-----#####-----#####
## D.Akinc - June 2017
## R code to visualize the Inverse Wishart, the Separation Strategy, the Scaled Inverse Wishart,
## and Huang Half-t priors
##
## Code is adapted from the code of Matt Simpson (2012)
## Url: http:// www.themattsimpson.com/2012/08/20/
##				prior-distributions-for-covariance-matrices-the-scaled-inverse-wishart-prior/ 
#####-----#####-----#####-----#####
#Required libraries:
require(MCMCpack)
require(ggplot2)
require(reshape)
require(gridExtra)
#install.packages(c("coda","mvtnorm","devtools","loo"))
#library(devtools)
#devtools::install_github("rmcelreath/rethinking")
library(rethinking)

#(1) Inverse Wishart Distribution
## simulates n samples from the IW prior
IW.test <- function(n, p, df, M){
  sig1 <- rep(0,n); sig2 <- rep(0,n); rho <- rep(0,n);
  for(i in 1:n){
    E <- riwish(df, M)
    sds <- sqrt(diag(E))
    rho[i] <- cov2cor(E)[2,1]; 
    sig1[i] <- sds[1]; sig2[i] <- sds[2];
  }
  return(data.frame(value=c(sig1, sig2, rho), 
                    parameter=c(rep("sigma1",n), rep("sigma2",n),rep("rho",n)), sam=rep(c(1:n),3)))
}

#(2.1) Separation Strategy with the setting of Barnard, et al. (2000)
# Correlation parameters from Inverse Wishart prior 
## simulates n samples from the SS prior
SS1.test <- function(n, p, b, s, df, M=diag(p)){
  sig1 <- rep(0,n)
  sig2 <- rep(0,n)
  rho <- rep(0,n)
  for(i in 1:n){
    E <- riwish(df, M)
    rho[i] <- cov2cor(E)[2,1]
    sds <- exp(mvrnorm(1,b, diag(s)))
    sig1[i] <- sds[1]
    sig2[i] <- sds[2]
  }
  return(data.frame(value=c(sig1, sig2, rho), parameter=c(rep("sigma1",n), rep("sigma2",n), rep("rho",n)), sam=rep(c(1:n),3)))
}

#(2.2) Separation Strategy with the setting of STAN manual (2016)
# Correlation parameters from LKJ prior
## simulates n samples from the SS prior
SS2.test <- function(n, p, bstan, sstan, eta){
  sig1 <- rep(0,n)
  sig2 <- rep(0,n)
  rho <- rep(0,n)
  for(i in 1:n){
    E <- rlkjcorr(1, p, eta)
    rho[i] <- cov2cor(E)[2,1]
    repeat {
      sds <- rcauchy(2,location=bstan,scale=sstan)
      if(sum(sds<=0)==0) {
        break
      }
    }
    sig1[i] <- sds[1]
    sig2[i] <- sds[2]
  }
  return(data.frame(value=c(sig1, sig2, rho), parameter=c(rep("sigma1",n), rep("sigma2",n), rep("rho",n)), sam=rep(c(1:n),3)))
}

#(3) Scaled Inverse Wishart Distribution 
## simulates a single sample from the sIW prior
sIW.sim <- function(p, b, s, df, M){
  R <- riwish(df, M)
  S <- diag(exp(mvrnorm(1,b, diag(s))))
  SRS <- S%*%R%*%S
  return(SRS)
}
## simulates n samples from the sIW prior
sIW.test <- function(n, p, b, s, df, M){
  sig1 <- rep(0,n); sig2 <- rep(0,n); rho <- rep(0,n);
  for(i in 1:n){
    E <- sIW.sim(2, b, s, df, M)
    sds <- sqrt(diag(E))
    rho[i] <- cov2cor(E)[2,1]; 
    sig1[i] <- sds[1]; sig2[i] <- sds[2]; 
  }
  return(data.frame(value=c(sig1, sig2, rho), 
                    parameter=c(rep("sigma1",n), rep("sigma2",n),rep("rho",n)), sam=rep(c(1:n),3)))
}

#(4) Huang Half-t Distribution
## simulates a single sample from the Hht prior
rhuangwand <- function(nu=2, a, A)
{
  if(missing(A)) k <- length(a)
  else {
    k <- length(A)
    a <- rgamma(k, 0.5, 1/A^2)}
  x <- riwish(nu + k - 1, 2*nu*diag(a))
  return(x)
}
## simulates n samples from the sIW prior
Hht.test <- function(n, p, nu=2, A=rep(1,p)){
  sig1 <- rep(0,n); sig2 <- rep(0,n); rho <- rep(0,n);
  for(i in 1:n){
    E <- rhuangwand(nu=2, A=A)
    sds <- sqrt(diag(E))
    rho[i] <- cov2cor(E)[2,1]; 
    sig1[i] <- sds[1]; sig2[i] <- sds[2];
    }
  return(data.frame(value=c(sig1, sig2, rho), 
                    parameter=c(rep("sigma1",n), rep("sigma2",n),rep("rho",n)), sam=rep(c(1:n),3)))
}

# Input values:
n <- 10000 #number of samples
p <- 2 #number of parameters
b <- rep(0,p) #b vector used in SS and SIW; mean of log-normal distribution for log of standard deviations
s <- rep(1,p) #xi vector used in SS and SIW; diag(s) is the variance of log-normal distribution for log of standard deviations
df <- p+1 #nu; degrees of freedom for IW, SS and Hht 
M <- diag(p) #p-dimensional scale matrix for all priors
nu <- 2 #nu; degrees of freedom for Hht
xi <- rep(1,p) #xi vector used in Hht
eta <- 1 #eta; shape parameter for SS2_STAN
bstan <- 1
sstan <- 1

#The prior distributions with the hyperparameters given in Table 1 and Alvarez, et al. (2014)
IWsam.1 <- IW.test(n, p, df, M)
SSsam.1 <- SS1.test(n, p, b=rep(log(0.72)/2,p), s, df, M)
SSsam.2 <- SS2.test(n, p, bstan=0.72*bstan, sstan=1*sstan, eta)
sIWsam.1 <- sIW.test(n, p, b, s, df, 0.8*M)
Hhtsam.1 <- Hht.test(n, p, nu, 1.04*xi)

IWsam.1$dens <- "1.IW"
SSsam.1$dens <- "2.SS"
SSsam.2$dens <- "3.SS_STAN"
sIWsam.1$dens <- "4.SIW"
Hhtsam.1$dens <- "5.Hht"


data.1 <- rbind(IWsam.1, SSsam.1, SSsam.2, sIWsam.1, Hhtsam.1)

#The plots given in the paper
par(mfrow=c(1,2))
#The Figure 3 (left panel): Plot of standard deviation
theme_set(theme_gray(base_size = 14))
pl3.1<-qplot(value, data=data.1[data.1$parameter=="sigma1",], bins=50, facets=dens~parameter, xlab="Standard Deviation",xlim=c(0,10))

#The Figure 3 (left panel):Plot of correlation coefficient
theme_set(theme_gray(base_size = 14))
pl3.2<-qplot(value, data=data.1[data.1$parameter=="rho",], facets=dens~parameter, xlab="Correlation")
grid.arrange(pl3.1, pl3.2, ncol=2)

data.1.melt <- melt(data.1, id=c("parameter", "dens", "sam"))
data.1.cast <- cast(data.1.melt, dens+sam~parameter)
# Figure 4: Plot of the first variance component vs. correlation coefficient
pl4<-qplot(sigma1^2, rho, data=data.1.cast, facets=.~dens, xlab="First variance component", ylab="Correlation", log="x")
pl4