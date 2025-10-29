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
  // real sigma;               // marginal std dev of GP
  real rho;        // lengthscale
  vector<lower=0>[C] lambda;     // horseshoe local
  array[C] real<lower=0> tau;        // observation noise sd
  // array[S,C] vector[N] gamma;    // subject specific observation noise
  // array[S,C] real<lower=-1,upper=1> rho_gamma;
  // real<lower=0> sigma_gamma;
  vector[N] f;                     // latent GP per channel
  real<lower=0> tau_sigma_rho;
  vector[C] alpha;
  // vector[C] delta;
  vector<lower=0>[C] beta;
  // vector<lower=0>[S] zeta;
  // real<lower=-1,upper=1> 
  vector<lower=-1, upper=1>[C] ar_phi;   // AR(1) coefficients per channel
}

transformed parameters {
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
  // vector[C] rho;
  // for (c in 1:C) rho[c] = -abs(phi[c]) * tau_sigma_rho * lambda[c] + mu_rho;
  
  vector<lower=0,upper=1>[C] phi = inv_logit(alpha);
  // vector[C] phi = alpha;
}

model {
  // Priors
  target += -0.5 * dot_self(alpha[node1] - alpha[node2]);
  sum(alpha) ~ normal(0, 0.01 * M);
  
  // target += -0.5 * dot_self(delta[node1] - delta[node2]);
  // sum(delta) ~ normal(0, 0.01 * C);
  // alpha ~ normal(0, tau_sigma_rho .* lambda);
  tau_sigma_rho ~ student_t(nut, 0, sigma_rho);
  lambda ~ student_t(nul, 0, phi);
  beta ~ student_t(nul, 0, tau_sigma_rho .* lambda);
  // lambda ~ cauchy(0, m);
  // sigma ~ normal(mu_sigma, sigma_sigma);
  // rho ~ normal(mu_rho, tau_sigma_rho .* lambda);
  // sigma ~ std_normal();
  rho ~ std_normal();
  // zeta ~ std_normal();
  // gamma[,1] ~ std_normal();
  // for ()
  // rho ~ normal(mu_rho, tau_sigma_rho .* sqrt(lambda2_tilde));
  // tau ~ normal(mu_tau, sigma_tau);
  // tau ~ std_normal();
  tau ~ normal(0, .05);
  ar_phi ~ normal(0, .1);
  // rho_gamma ~ normal(0, .5);
  // sigma_gamma ~ std_normal();
  // gamma[,1] ~ std_normal();
  // for (c in 1:C) {
  //   for (s in 1:S) {
  //     rho_gamma[s,c] ~ normal(0, .5);
  //     gamma[s,c,1] ~ std_normal();
  //     gamma[s,c,2:N] ~ normal(rho_gamma[s,c] * gamma[s,c,1:(N-1)], sigma_gamma);
  //   }
  // }
  // alpha ~ normal(0, m);
  
  // for (c in 1:C) {
    target += gp_graph_exp_quad_cov_lpdf(f | zeros_vector(N), x, 
                                         1, 
                                         // exp(mu_rho + 
                                         // tau_sigma_rho * lambda[c] * rho[c]), 
                                         exp(rho),
                                         edge_index);
  // }

  // Likelihood: vectorized over subjects
  
  // for (c in 1:C) {
  //   matrix[S, N] y_c = to_matrix(y[, , c]);
  //   to_vector(y_c) ~ normal(rep_vector(beta[c], S * N) .* rep_vector(f, S), tau[c]);
  // }
  
  for (c in 1:C) {
    for (s in 1:S) {
      vector[N] mu = beta[c] * f;
      vector[N] resid = to_vector(y[s, , c]) - mu;

    // AR(1) model on residuals
      target += normal_lpdf(resid[1] | 0, tau[c] / sqrt(1 - square(ar_phi[c])));
      for (t in 2:N) {
        target += normal_lpdf(resid[t] | ar_phi[c] * resid[t - 1], tau[c]);
      }
    }
  }
  
  // for (c in 1:C) {
  //   for (s in 1:S) {
  //     for (n in 1:N) {
  //     // elementwise independent normal
  //       y[s,n,c] ~ normal(beta[c] * f[n], tau[c]);
  //     }
  //   }
  // }
}


generated quantities {
  array[C,S] vector[N] mean_res;         // your existing predictions
  array[C,S] vector[N] resid_out;     // store residuals for each subject & channel
  array[C] vector[N] pred_res;
  array[C,S] vector[N] preds_res;
  // for (s in 1:S) {
  for (c in 1:C) {
    pred_res[c] = beta[c] * f;

    for (s in 1:S) {
      mean_res[c,s] = beta[c] * f; 
      vector[N] mu = beta[c] * f;
      vector[N] resid = to_vector(y[s, , c]) - mu;
      resid_out[c, s] = resid;         // save residuals
      // pred_res[c, s] = pred_res[c] + resid_out[c, s];
    }
  }
  // }
  
  // array[C, S] vector[N] yhat;   // overall predicted signal per channel and subject
  // 
  for (c in 1:C) {
    for (s in 1:S) {
      preds_res[c,s,1] = mean_res[c,s,1] + normal_rng(0, 
                                                tau[c] / 
                                                  sqrt(1 - square(ar_phi[c])));
      for (n in 2:N) {
        preds_res[c,s,n] = mean_res[c,s,n] + 
                              normal_rng(ar_phi[c]*resid_out[c,s,n-1], 
                                         tau[c]);
      }
    }
  }
}



















