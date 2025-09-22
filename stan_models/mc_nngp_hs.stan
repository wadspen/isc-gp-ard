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
  real r;                            // horsehoe regularization
}

parameters {
  real sigma;               // marginal std dev of GP
  vector[C] rho;        // lengthscale
  vector<lower=0>[C] lambda;     // horseshoe local
  array[C] real<lower=0> tau;        // observation noise sd
  matrix[N, C] f;                     // latent GP per channel
}

// transformed parameters {
//   // matrix[N, C] f;                     // latent GP per channel
//   //
//   // for (c in 1:C) {
//   //   // Compute GP once per channel
//   //   f[,c] = gp_graph_exp_quad_cov(zeros_vector(N), x, sigma, 1/rho[c], edge_index);
//   // }
//   vector<lower=0>[C] lambda2_tilde;
// 
//   for (c in 1:C) lambda2_tilde[c] = (r^2 * lambda[c]^2) /
//                                     (r^2 + sigma_rho^2 * lambda[c]^2);
// }

model {
  // Priors
  lambda ~ normal(0, 1);
  sigma ~ normal(mu_sigma, sigma_sigma);
  rho   ~ normal(mu_rho, sigma_rho .* lambda);
  // rho ~ normal(mu_rho, sqrt(lambda2_tilde));
  tau   ~ normal(mu_tau, sigma_tau);
  
  for (c in 1:C) {
    target += gp_graph_matern52_cov_lpdf(f[, c] | zeros_vector(N), x, 
                                         exp(sigma), 
                                         exp(rho[c]), edge_index);
  }

  // Likelihood: vectorized over subjects
  for (c in 1:C) {
    for (s in 1:S) {
      // elementwise independent normal
      y[s, , c] ~ normal(f[, c], tau[c]);
    }
  }
}
