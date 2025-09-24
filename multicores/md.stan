data{
  int N;
  vector[N] x;
}

parameters{
  real meanV;
}

model{
  // priors
  meanV ~ normal(0, 1);
  
  // model
  vector[N] sub;
  
  for (n in 1:N){
    sub[n] = 1 + x[n]/2;
  }
  sub ~ normal(meanV, 2);
}
