functions {
  #include gptools/util.stan
  #include gptools/fft.stan
  #include gptools/graph.stan
}

data {
  int<lower=1> N;                    // number of grid points
  int<lower=1> S;                    // number of subjects
  int<lower=1> C;                    // number of channels / parcels
  array[S, N, C] real y;             // observations
  array[N] vector[1] x;              // grid locations
  array[2, N - 1] int edge_index;    // graph edges for GP
  real mu_rho;
  real<lower=0> sigma_rho;
  real mu_sigma;
  real<lower=0> sigma_sigma;
  real mu_tau;
  real<lower=0> sigma_tau;
}

parameters {
  real<lower=0> sigma;               // marginal std dev of GP
  array[C] real<lower=0> rho;        // lengthscale
  array[C] real<lower=0> tau;        // observation noise sd
  matrix[N, C] f;                     // latent GP per channel
}

// transformed parameters {
//   matrix[N, C] f;                     // latent GP per channel
// 
//   for (c in 1:C) {
//     // Compute GP once per channel
//     f[,c] = gp_graph_exp_quad_cov(zeros_vector(N), x, sigma, 1/rho[c], edge_index);
//   }
// }

model {
  // Priors
  sigma ~ normal(mu_sigma, sigma_sigma);
  rho   ~ normal(mu_rho, sigma_rho);
  tau   ~ normal(mu_tau, sigma_tau);
  
  for (c in 1:C) {
    target += gp_graph_exp_quad_cov_lpdf(f[, c] | zeros_vector(N), x, sigma, 1/rho[c], edge_index);
  }

  // Likelihood: vectorized over subjects
  for (c in 1:C) {
    for (s in 1:S) {
      // elementwise independent normal
      y[s, , c] ~ normal(f[, c], tau[c]);
    }
  }
}
