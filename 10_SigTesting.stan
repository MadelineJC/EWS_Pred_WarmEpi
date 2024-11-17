data {
  int <lower = 0> N;
  int <lower = 0> y[N];
}

parameters {
  vector[N] mu; // Mean value
  real <lower = 0> sigma; // Variation
  real ac; // Autocorrelation
}

model {
  mu[1] ~ lognormal(0, 1); // Prior on mu
  for (i in 2:N){
    mu[i] ~ normal(ac * mu[i - 1], sigma); // Autoregressive part
  }
  y ~ poisson_log(mu); // Vectorized over the N elements
}

generated quantities {
  int <lower = 0> y_tilde[N];
  for (i in 1:N){
    y_tilde[i] = poisson_log_rng(mu[i]);
  }
}
