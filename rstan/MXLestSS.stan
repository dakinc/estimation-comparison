//#####-----#####-----#####-----#####
//## D.Akinc - June 2017
//##
//## (2) STAN code for mixed logit estimation with the Separation Stategy 
//## prior of Barnard et al.(2000) for the covariance matrix
//## File: 'MXLestSS.stan'
//##
//## The code is adapted from Guhl D. & Gabel S. (2016). Discrete Choice
//## Models in R[Presentation]. Retrieved from 
//## Url: https:// github.com/eRum2016/Presentations-participants/blob/master/
//## 			13.10/Methodology%202/eRum_GG2016.pdf
//#####-----#####-----#####-----#####

data {
  int<lower=1> S;                           // number of choice sets
  int<lower=1> P;                           // number of parameters (utility)
  int<lower=1> D;                           // number of parameters (beta)
  int<lower=2> K;                           // number of alternatives
  int<lower=1> N;                           // number of respondents
  int<lower=1,upper=K> y[S];                // alternative selected
  int<lower=1,upper=N> id[S];               // respondent index
  matrix[K,P] x[S];                         // array of design matrices
  matrix[N,D] d;                            // matrix of demographics
}

parameters {
  vector[P] mu;                             // means for beta
  corr_matrix[P] Omega;                		// cholesky of correlation matrix
  vector<lower=0>[P] sig;                   // variances
  vector[P] beta[N];                        // individual level parameters

}

transformed parameters {
  cov_matrix[P] Sigma;                      // cholesky factors
  vector[K] v[S];                           // utility vector for each choice set
  matrix[P,P] Vprior;                       // V for Inverse Wishart

  Sigma = quad_form_diag(Omega, sig);

  for (s in 1:S) {
    v[s] = x[s] * beta[id[s]];
  }
  Vprior = (P+1) *diag_matrix(rep_vector(1, P));
}

model {
  mu ~ normal(0, 100);                 		// prior for mean parameters
  sig ~ lognormal(0, 1);                    // prior for sddev
  Omega ~ inv_wishart((P+1),Vprior);        // prior for corr - flat prior

  for (n in 1:N){
    beta[n] ~ multi_normal(mu,Sigma);
  }

  for (s in 1:S) {
    y[s] ~ categorical_logit(v[s]);
  }
}

generated quantities {
  real<upper=0> log_lik;                    // log_lik for each sample

  log_lik = 0;
  for (s in 1:S) {
    log_lik = log_lik + categorical_logit_lpmf(y[s] | v[s]);
  }
}
