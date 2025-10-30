//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;                           //num observations
  int<lower=0> S;                           //num subjects
  int<lower=0> R;                           //num ROIs/voxels
  int<lower=0> V;                           //num obs per ROI/voxel
  array[V] int n1_inds;                        //indices for xi n1
  array[V] int n2_inds;                        //indices for xi n2
  matrix[V,R] y;
  vector[S] subjects;
  array[N] int n1s;
  array[N] int n2s;
  array[R] int ROIs;
  array[N] int rs;
  matrix[V,V] W;
  matrix[V,V] I;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real alpha0;
  real<lower=0> sigma;
  vector[S] xi;
  vector[R] pi0;
  real<lower=0> lambda;
  real<lower=0> tau;
}

transformed parameters {
  real rho = (lambda^2 + tau^2) / (2 * lambda^2 + tau^2 + sigma^2);
  matrix[V,V] cov = square(sigma) * (rho * W + I) + 1e-6*I;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  lambda ~ std_normal();
  tau ~ std_normal();
  xi ~ normal(0, lambda);
  pi0 ~ normal(0, tau);
  
  
  for (r in 1:R) {
    y[,r] ~ multi_normal(alpha0 + xi[n1_inds] + xi[n2_inds] + pi0[r], 
                            cov);
  }
}












