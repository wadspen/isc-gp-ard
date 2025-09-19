functions {
  #include gptools/util.stan
  #include gptools/fft.stan
  #include gptools/graph.stan
}

data {
  int<lower=1> N;                    // number of regularly-spaced grid points
  int<lower=1> S;                    // number of subjects
  int<lower=1> C;
  // vector[C] coords;
  // vector[n] y;                       // observations on the grid
  array[S, N, C] real y;
  array[N] vector[1] x;
  array[2, N - 1] int edge_index;
  real mu_rho;
  real<lower=0> sigma_rho;
  real mu_sigma;
  real<lower=0> sigma_sigma;
  real mu_tau;
  real<lower=0> sigma_tau;
  // real<lower=0> sigma;
}

parameters {
  real<lower=0> sigma;               // marginal std dev of GP
  array[C] real<lower=0> rho;                 // lengthscale (in same units as dx)
  array[C] real<lower=0> tau;                 // observation noise sd
  // vector[N] z;
  matrix[N, C] f;
}

// transformed parameters {
//   vector[N] f = gp_graph_exp_quad_cov(z, zeros_vector(N), x, sigma, 
//                                           rho, edge_index);
//   
// }

model {
  // priors (tune these to the problem)
  sigma ~ normal(mu_sigma, sigma_sigma);
  rho   ~ normal(mu_rho, sigma_rho);    // or lognormal, or half-normal centered at a reasonable scale
  tau   ~ normal(mu_tau, sigma_tau);
  // z ~ std_normal();
  
  
  // likelihood
  for (c in 1:C) {
    f[,c] ~ gp_graph_exp_quad_cov(zeros_vector(N), x, sigma, 
                                          1/rho[c], edge_index);
    for (s in 1:S) {
      to_vector(y[s,,c]) ~ multi_normal(f[,c], 
                             diag_matrix(rep_vector(square(tau[c]), N)));
    }
  }
}
