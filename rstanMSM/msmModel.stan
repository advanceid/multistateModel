// Build a multistate model with rstan to estimate transition intensity
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
  real<lower=0> a12;
  real<lower=0> a13;
  real<lower=0> a15;
  
  real<lower=0> a21;
  real<lower=0> a24;
  real<lower=0> a25;
  
  real<lower=0> a31;
  real<lower=0> a34;
  real<lower=0> a35;
  
  real<lower=0> a42;
  real<lower=0> a43;
  real<lower=0> a45;
}

transformed parameters {
  matrix[5, 5] mat;
  
  mat[1,1] = -a12 - a13 - a15;
  mat[1,2] = a12;
  mat[1,3] = a13;
  mat[1,4] = 0;
  mat[1,5] = a15;
  
  mat[2,1] = a21;
  mat[2,2] = -a21 - a24 - a25;
  mat[2,3] = 0;
  mat[2,4] = a24;
  mat[2,5] = a25;
  
  mat[3,1] = a31;
  mat[3,2] = 0;
  mat[3,3] = -a31 - a34 - a35;
  mat[3,4] = a34;
  mat[3,5] = a35;
  
  mat[4,1] = 0;
  mat[4,2] = a42;
  mat[4,3] = a43;
  mat[4,4] = -a42 - a43 - a45;
  mat[4,5] = a45;
  
  mat[5,1] = 0;
  mat[5,2] = 0;
  mat[5,3] = 0;
  mat[5,4] = 0;
  mat[5,5] = 0;
}

model {
  // Priors
  a12 ~ exponential(1);
  a13 ~ exponential(1);
  a15 ~ exponential(1);
  
  a21 ~ exponential(1);
  a24 ~ exponential(1);
  a25 ~ exponential(1);
  
  a31 ~ exponential(1);
  a34 ~ exponential(1);
  a35 ~ exponential(1);
  
  a42 ~ exponential(1);
  a43 ~ exponential(1);
  a45 ~ exponential(1);
  
  // Likelihood
  for (n in 1:N) {
    matrix[5, 5] prob_n = matrix_exp(mat * (t1[n] - t0[n]));
    s1[n] ~ categorical(prob_n[s0[n],]');
  }
}
