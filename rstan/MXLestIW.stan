//#####-----#####-----#####-----#####
//## D.Akinc - June 2017
//##
//## (1) STAN code for mixed logit estimation with the Inverse Wishart prior 
//## for the covariance matrix
//## File: 'MXLestIW.stan
//##
//## The code is adapted from Guhl D. & Gabel S. (2016). Discrete Choice
//## Models in R[Presentation]. Retrieved from 
//## Url: https:// github.com/eRum2016/Presentations-participants/blob/master/
//## 			13.10/Methodology%202/eRum_GG2016.pdf
//#####-----#####-----#####-----#####

data {
  int<lower=1> S;                    // number of choice sets
  int<lower=1> P;                    // number of parameters (utility)
  int<lower=1> D;                    // number of parameters (beta)
  int<lower=2> K;                    // number of alternatives
  int<lower=1> N;                    // number of respondents
  int<lower=1,upper=K> y[S];         // choice data
  int<lower=1,upper=N> id[S];        // respondent index
  matrix[K,P] x[S];                  // Design matrices
  matrix[N,D] d;                     // matrix of demographics
}

parameters {
  vector[P] mu;                      // MU_BETA
  cov_matrix[P] Sigma;               // SIGMA_BETA
  vector[P] beta[N];                 // BETA_N
}

transformed parameters {
  vector[K] v[S];                    // utility vector for each choice set
  matrix[P,P] Vprior;                // scale parameter for inverse Wishart  
  for (s in 1:S) {
    v[s] = x[s] * beta[id[s]];
	}
  Vprior = (P+1) * diag_matrix(rep_vector(1, P));
}

model {
  mu ~ normal(0, 100);                // prior for MU_BETA
  Sigma ~ inv_wishart((P+1),Vprior);  // prior for SIGMA_BETA
  for (n in 1:N){
    beta[n] ~ multi_normal(mu,Sigma);
  }
  for (s in 1:S) {
    y[s] ~ categorical_logit(v[s]);
  }
}

generated quantities {
  real<upper=0> log_lik;              // log_likelihood for each sample
  log_lik = 0;
  for (s in 1:S) {
    log_lik = log_lik + categorical_logit_lpmf(y[s] | v[s]);
  }
}
