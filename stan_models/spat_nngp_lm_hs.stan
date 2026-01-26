functions {
  #include gptools/util.stan
  #include gptools/fft.stan
  #include gptools/graph.stan
}

data {
  int<lower=1> N;                    // number of grid points time
  int<lower=1> M;                    // number of grid points space
  int<lower=0> n_edges;              // number of edges on spatial domain
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
  real<lower=0> r;                            // horsehoe regularization
  real<lower=0> nut;
  real<lower=0> nul;
  real<lower=0> m;
  array[n_edges] int<lower=1, upper=M> node1; // node1[i] adjacent to node2[i]
  array[n_edges] int<lower=1, upper=M> node2; // and node1[i] < node2[i]
}

parameters {
  real rho;        // lengthscale
  vector<lower=0>[C] lambda;     // horseshoe local
  array[C] real<lower=0> tau;        // observation noise sd
  vector[N] f;                     // latent GP per channel
  real<lower=0> tau_sigma_rho;
  vector[C] alpha;
  vector<lower=0>[C] beta;
}

transformed parameters {
  vector<lower=0,upper=1>[C] phi = inv_logit(alpha);
}

model {
  // Priors
  target += -0.5 * dot_self(alpha[node1] - alpha[node2]);
  sum(alpha) ~ normal(0, 0.01 * M);
  tau_sigma_rho ~ student_t(nut, 0, sigma_rho);
  lambda ~ student_t(nul, 0, phi);
  beta ~ student_t(nul, 0, tau_sigma_rho .* lambda);
  rho ~ std_normal();
  tau ~ normal(mu_tau, sigma_tau);
    target += gp_graph_exp_quad_cov_lpdf(f | zeros_vector(N), x, 
                                          1, 
                                         exp(rho),
                                         edge_index);

  
  for (c in 1:C) {
    for (s in 1:S) {
      y[s, , c] ~ normal(beta[c] * f, tau[c]);
    }
  }
  

}

// generated quantities {
//   array[C] vector[N] pred_res;
//   for (c in 1:C) {
//     pred_res[c] = beta[c]*f;
//   }
//   
//   
//   
// }



















