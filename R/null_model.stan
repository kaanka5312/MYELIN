data {
  int<lower=1> N;              // Total number of observations
  vector[N] y;                 // Response variable for all observations
}

transformed data{
  real y_mean;                                   // mean of y; see mu prior
  
  y_mean = mean(y); 
}

parameters {
  real mu;                     // Common mean for all observations
  real<lower=0> sigma;         // Common standard deviation for all observations
  real<lower=0, upper=100> nu; // Degrees of freedom for the Student's t-distribution
}

model {
  // Priors
  mu ~ normal(y_mean, 10);          // Assuming a wide prior centered at 0 for the mean
  sigma ~ cauchy(0, 5);        // Using a Cauchy prior for the standard deviation
  nu    ~ exponential(1.0/29); // A weakly informative prior for degrees of freedom
  
  // Likelihood
  y ~ student_t(nu, mu, sigma); // Assuming all data points come from a single Student's t-distribution
}

generated quantities {
  vector[N] log_lik;           // log likelihood for each observation
  
  for (n in 1:N) {
    log_lik[n] = student_t_lpdf(y[n] | nu, mu, sigma);
  }
}


