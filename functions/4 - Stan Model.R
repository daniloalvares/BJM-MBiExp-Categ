# JOINT MODEL OF MULTIPLE LONGITUDINAL AND CATEGORICAL OUTCOMES

JM_Biexp_Cat_FIT <- write_stan_file("
functions {
  // Bi-exponential specification
  vector nonlinear_predictor(array[] int ID, vector times, vector theta, matrix bi){
     int N = num_elements(times);
     vector[N] logBi = theta[1] + bi[ID,1];      // log-baseline
     vector[N] Gi = exp(theta[2] + bi[ID,2]);    // growth rate
     vector[N] Di = exp(theta[3] + bi[ID,3]);    // decay rate

     vector[N] out = logBi + log( exp(rows_dot_product(Gi, times)) + exp(-rows_dot_product(Di, times)) - 1 );

     return out;
  }
}
data {
  int N_M;
  int N_F;
  int n_M;
  int n_F;
  int n;
  int nbeta;
  array[N_M] real y_M;
  array[N_F] real y_F;
  array[N_M] int<lower=1,upper=n_M> ID_M1;
  array[N_M] int<lower=1,upper=n> ID_M2;
  array[N_F] int<lower=1,upper=n_F> ID_F1;
  array[N_F] int<lower=1,upper=n> ID_F2;
  vector[N_M] times_M;
  vector[N_F] times_F;
  array[n] int z;
  matrix[n,nbeta] X;
}
transformed data {
  vector[nbeta] zerosb = rep_vector(0, nbeta);
  vector[6] zerosa = rep_vector(0, 6);
}
parameters {
  vector[3] theta_M;
  vector[3] theta_F;
  real<lower=0> sigma2_M;
  real<lower=0> sigma2_F;
  cov_matrix[3] Omega_M;
  cov_matrix[3] Omega_F;
  matrix[nbeta,2] beta_raw;
  matrix[6,2] alpha_raw;
  matrix[n_M,3] bi_M;
  matrix[n_F,3] bi_F;
}
transformed parameters {
  // log-baseline, log-growth and log-decay paramaters
  matrix[n,2] lBi;
  matrix[n,2] lGi;
  matrix[n,2] lDi;
  // Coefficients for the categorical model
  matrix[nbeta,3] beta = append_col(beta_raw, zerosb);
  matrix[6,3] alpha = append_col(alpha_raw, zerosa);
  
  // M-spike
  lBi[,1] = rep_vector(theta_M[1],n);
  lBi[ID_M2,1] = theta_M[1] + bi_M[ID_M1,1];
  lGi[,1] = rep_vector(theta_M[2],n);
  lGi[ID_M2,1] = theta_M[2] + bi_M[ID_M1,2];
  lDi[,1] = rep_vector(theta_M[3],n);
  lDi[ID_M2,1] = theta_M[3] + bi_M[ID_M1,3];
  
  // FLC
  lBi[,2] = rep_vector(theta_F[1],n);
  lBi[ID_F2,2] = theta_F[1] + bi_F[ID_F1,1];
  lGi[,2] = rep_vector(theta_F[2],n);
  lGi[ID_F2,2] = theta_F[2] + bi_F[ID_F1,2];
  lDi[,2] = rep_vector(theta_F[3],n);
  lDi[ID_F2,2] = theta_F[3] + bi_F[ID_F1,3];
}
model {
  // BI-EXPONENTIAL SPECIFICATION
  vector[N_M] nonlinpred_M = nonlinear_predictor(ID_M1, times_M, theta_M, bi_M);
  vector[N_F] nonlinpred_F = nonlinear_predictor(ID_F1, times_F, theta_F, bi_F);
  // CATEGORICAL SPECIFICATION
  matrix[n,3] linpred = X * beta + append_col(append_col(lBi, lGi), lDi) * alpha;
   
  // LONGITUDINAL NORMAL LOG-LIKELIHOOD
  target += normal_lpdf(y_M | nonlinpred_M, sqrt(sigma2_M));
  target += normal_lpdf(y_F | nonlinpred_F, sqrt(sigma2_F));

  // CATEGORICAL LOG-LIKELIHOOD
  for(i in 1:n){ target += categorical_logit_lpmf(z[i] | to_vector(linpred[i,])); }
  
  // LOG-PRIORS
  // Longitudinal fixed effects
  target += normal_lpdf(theta_M | 0, 10);
  target += normal_lpdf(theta_F | 0, 10);
   
  // Residual error variance
  target += cauchy_lpdf(sigma2_M | 0, 5);
  target += cauchy_lpdf(sigma2_F | 0, 5);
  
  // Random-effects variance-covariance matrices
  target += inv_wishart_lpdf(Omega_M | 4, diag_matrix(rep_vector(1,3)));
  target += inv_wishart_lpdf(Omega_F | 4, diag_matrix(rep_vector(1,3)));  
  
  // Random-effects
  for(i in 1:n_M){ target += multi_normal_lpdf(bi_M[i,1:3] | rep_vector(0,3), Omega_M); }
  for(i in 1:n_F){ target += multi_normal_lpdf(bi_F[i,1:3] | rep_vector(0,3), Omega_F); }  
  
  target += normal_lpdf(to_vector(beta) | 0, 10);
  target += normal_lpdf(to_vector(alpha) | 0, 10);

}")


JM_Biexp_Cat_VI <- write_stan_file("
functions {
  // Bi-exponential specification
  vector nonlinear_predictor(array[] int ID, vector times, vector theta, matrix bi){
     int N = num_elements(times);
     vector[N] logBi = theta[1] + bi[ID,1];      // log-baseline
     vector[N] Gi = exp(theta[2] + bi[ID,2]);    // growth rate
     vector[N] Di = exp(theta[3] + bi[ID,3]);    // decay rate

     vector[N] out = logBi + log( exp(rows_dot_product(Gi, times)) + exp(-rows_dot_product(Di, times)) - 1 );

     return out;
  }
}
data {
  int N_M;
  int N_F;
  int n_M;
  int n_F;
  int n;
  int nbeta;
  array[N_M] real y_M;
  array[N_F] real y_F;
  array[N_M] int<lower=1,upper=n_M> ID_M1;
  array[N_M] int<lower=1,upper=n> ID_M2;
  array[N_F] int<lower=1,upper=n_F> ID_F1;
  array[N_F] int<lower=1,upper=n> ID_F2;
  vector[N_M] times_M;
  vector[N_F] times_F;
  array[n] int z;
  matrix[n,nbeta] X;
  array[n] int<lower=1,upper=n> IMP_B_M;
  array[n] int<lower=1,upper=n> IMP_B_F;
  array[n] int<lower=1,upper=n> IMP_G_M;
  array[n] int<lower=1,upper=n> IMP_G_F;
  array[n] int<lower=1,upper=n> IMP_D_M;
  array[n] int<lower=1,upper=n> IMP_D_F;
}
transformed data {
  vector[nbeta] zerosb = rep_vector(0, nbeta);
  vector[6] zerosa = rep_vector(0, 6);
}
parameters {
  vector[3] theta_M;
  vector[3] theta_F;
  real<lower=0> sigma2_M;
  real<lower=0> sigma2_F;
  cov_matrix[3] Omega_M;
  cov_matrix[3] Omega_F;
  matrix[nbeta,2] beta_raw;
  matrix[6,2] alpha_raw;
  matrix[n_M,3] bi_M;
  matrix[n_F,3] bi_F;
}
transformed parameters {
  // log-baseline, log-growth and log-decay paramaters
  matrix[n,2] lBi;
  matrix[n,2] lGi;
  matrix[n,2] lDi;
  // Coefficients for the categorical model
  matrix[nbeta,3] beta = append_col(beta_raw, zerosb);
  matrix[6,3] alpha = append_col(alpha_raw, zerosa);
  
  // M-spike
  lBi[,1] = rep_vector(theta_M[1],n);
  lBi[ID_M2,1] = theta_M[1] + bi_M[ID_M1,1];
  lGi[,1] = rep_vector(theta_M[2],n);
  lGi[ID_M2,1] = theta_M[2] + bi_M[ID_M1,2];
  lDi[,1] = rep_vector(theta_M[3],n);
  lDi[ID_M2,1] = theta_M[3] + bi_M[ID_M1,3];
  
  // FLC
  lBi[,2] = rep_vector(theta_F[1],n);
  lBi[ID_F2,2] = theta_F[1] + bi_F[ID_F1,1];
  lGi[,2] = rep_vector(theta_F[2],n);
  lGi[ID_F2,2] = theta_F[2] + bi_F[ID_F1,2];
  lDi[,2] = rep_vector(theta_F[3],n);
  lDi[ID_F2,2] = theta_F[3] + bi_F[ID_F1,3];
}
generated quantities {
  vector[N_M] nonlinpred_M = nonlinear_predictor(ID_M1, times_M, theta_M, bi_M);
  vector[N_F] nonlinpred_F = nonlinear_predictor(ID_F1, times_F, theta_F, bi_F);
  matrix[n,3] linpred = X * beta + 
                        append_col(append_col(append_col(lBi[IMP_B_M,1], lBi[IMP_B_F,2]), 
                                              append_col(lGi[IMP_G_M,1], lGi[IMP_G_F,2])),
                                              append_col(lDi[IMP_D_M,1], lDi[IMP_D_F,2])) * alpha;
  matrix[n,2] lB = lBi;
  matrix[n,2] lG = lGi;
  matrix[n,2] lD = lDi;
  matrix[n,3] probs;
  vector[n] log_lik;
  for(i in 1:n){
     probs[i,] = softmax(to_vector(linpred[i,]))';
     log_lik[i] = categorical_logit_lpmf(z[i] | to_vector(linpred[i,]));
  }
}")
