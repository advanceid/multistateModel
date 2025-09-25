// Build a multistate model with rstan to estimate transition intensity
// Estimating log transition intensity instead of original form
data{
  int<lower = 1> N;   // number of subjects
  int<lower = 1> P;   // number of predictors
  
  matrix[N, P]   X;   // predictors
  vector[N]      t0;   // starting time
  vector[N]      t1;   // ending time
  
  int<lower=1, upper=5> s0[N];   // starting state
  int<lower=1, upper=5> s1[N];   // ending state
  
}


parameters{
  // real<lower=0> a12;
  // real<lower=0> a13;
  // real<lower=0> a15;
  // 
  // real<lower=0> a21;
  // real<lower=0> a24;
  // real<lower=0> a25;
  // 
  // real<lower=0> a31;
  // real<lower=0> a34;
  // real<lower=0> a35;
  // 
  // real<lower=0> a42;
  // real<lower=0> a43;
  // real<lower=0> a45;
  
  real loga12;
  real loga13;
  real loga15;
  
  real loga21;
  real loga24;
  real loga25;

  real loga31;
  real loga34;
  real loga35;

  real loga42;
  real loga43;
  real loga45;
}

transformed parameters {
  matrix[5, 5] mat;
  
  mat[1,1] = -exp(loga12) - exp(loga13) - exp(loga15);
  mat[1,2] = exp(loga12);
  mat[1,3] = exp(loga13);
  mat[1,4] = 0;
  mat[1,5] = exp(loga15);
  
  mat[2,1] = exp(loga21);
  mat[2,2] = -exp(loga21) - exp(loga24) - exp(loga25);
  mat[2,3] = 0;
  mat[2,4] = exp(loga24);
  mat[2,5] = exp(loga25);
  
  mat[3,1] = exp(loga31);
  mat[3,2] = 0;
  mat[3,3] = -exp(loga31) - exp(loga34) - exp(loga35);
  mat[3,4] = exp(loga34);
  mat[3,5] = exp(loga35);
  
  mat[4,1] = 0;
  mat[4,2] = exp(loga42);
  mat[4,3] = exp(loga43);
  mat[4,4] = -exp(loga42) - exp(loga43) - exp(loga45);
  mat[4,5] = exp(loga45);
  
  mat[5,1] = 0;
  mat[5,2] = 0;
  mat[5,3] = 0;
  mat[5,4] = 0;
  mat[5,5] = 0;
}

model {
  // Priors
  loga12 ~ normal(0, 1);
  loga13 ~ normal(0, 1);
  loga15 ~ normal(0, 1);
  
  loga21 ~ normal(0, 1);
  loga24 ~ normal(0, 1);
  loga25 ~ normal(0, 1);
  
  loga31 ~ normal(0, 1);
  loga34 ~ normal(0, 1);
  loga35 ~ normal(0, 1);
  
  loga42 ~ normal(0, 1);
  loga43 ~ normal(0, 1);
  loga45 ~ normal(0, 1);
  
  // Likelihood
  for (n in 1:N) {
    matrix[5, 5] prob_n = matrix_exp(mat * (t1[n] - t0[n]));
    s1[n] ~ categorical(prob_n[s0[n],]');
  }
}
