functions {
  #include gptools/util.stan
  #include gptools/fft.stan
}

data {
  int<lower=1> N;                    // number of regularly-spaced grid points
  int<lower=1> S;                    // number of subjects
  // vector[n] y;                       // observations on the grid
  matrix[S, N] y;
  real<lower=0> dx;                  // grid spacing (in x-units)
  int<lower=1> nf;                   // length of RFFT = n%/%2 + 1
  // vector[nf] omega;                  // angular frequencies (omega_k = 2*pi * freq_k) for k=0..n/2
  real mu_rho;
  real<lower=0> sigma_rho;
  real mu_sigma;
  real<lower=0> sigma_sigma;
  real mu_tau;
  real<lower=0> sigma_tau;
  real<lower=0> sigma;
}

parameters {
  // real<lower=0> sigma;               // marginal std dev of GP
  real<lower=0> rho;                 // lengthscale (in same units as dx)
  real<lower=0> tau;                 // observation noise sd
  // vector[n] zf;                     // white noise in Fourier domain (non-centered)
  vector[N] f;
}

transformed parameters {
  // cov_rfft is the real-FFT of the covariance kernel evaluated on the domain
  vector[nf] cov_rfft;

  // --- squared-exponential spectral density (continuous form) ---
  // For an SE kernel k(t)=sigma^2 * exp(- t^2 / (2 rho^2)),
  // the continuous spectral density (as function of angular freq omega) is
  //   S(omega) = sigma^2 * sqrt(2*pi) * rho * exp(-0.5 * (rho * omega)^2)
  //
  // We use that spectral shape sampled at the discrete angular frequencies omega[k].
  // (If you use a different convention for omega, adapt accordingly.)
  // for (k in 1:nf) {
  //   cov_rfft[k] = sigma^2 * sqrt(2 * pi()) * rho * exp(-0.5 * square(rho * omega[k]));
  //   // small jitter for numeric stability:
  //   cov_rfft[k] = cov_rfft[k] + 1e-10;
  // }
  cov_rfft = gp_periodic_exp_quad_cov_rfft(N, sigma, 1/rho, N)
        + 1e-9;
}

model {
  // priors (tune these to the problem)
  // sigma ~ normal(mu_sigma, sigma_sigma);
  rho   ~ normal(mu_rho, sigma_rho);    // or lognormal, or half-normal centered at a reasonable scale
  tau   ~ normal(mu_tau, sigma_tau);

  // Non-centered parameterization: zf ~ normal(0, 1) makes gp coefficients white.
  // zf ~ normal(0, 1);

  // Transform zf to f (real-space GP) via gp_rfft_inv or using gp_rfft helpers.
  // gptools provides gp_rfft and helpers that transform between f and its Fourier whitened representation.
  // Here we call gp_rfft_inv-like functionality by using gp_rfft_log_abs_det_jac and gp_rfft.
  // The typical pattern: f = gp_rfft_inv(zf, cov_rfft). gptools exposes gp_rfft transform functions.
  f ~ gp_rfft(zeros_vector(N), cov_rfft); // NOTE: function name depends on gptools version; see docs.

  // likelihood
  for (s in 1:S) {
    y[s,] ~ normal(f, tau);
  }
}
