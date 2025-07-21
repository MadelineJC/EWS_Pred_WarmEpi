data {
  int <lower = 0> N;
  int <lower = 0> y[N];
}

parameters {
  vector[N] mu; // Mean value
  real <lower = 0> sigma; // Variation
  real ac; // Autocorrelation
}

transformed parameters {
  vector[N] alpha;
  vector[N] beta;
  for (i in 1:N){
    alpha[i] = log(mu[i])^2 / sigma^2;
    beta[i] = log(mu[i]) / sigma^2;
  }
}

model {
  mu[1] ~ lognormal(0, 1); // Prior on mu
  for (i in 2:N){
    mu[i] ~ normal(ac * mu[i - 1], sigma); // Autoregressive part
  }
  y ~ neg_binomial(alpha, beta); // Vectorized over the N elements
}

generated quantities {
  int <lower = 0> y_tilde[N];
  for (i in 1:N){
    y_tilde[i] = neg_binomial_rng(alpha[i], beta[i]);
  }
}
