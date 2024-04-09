data {
  int<lower = 1> N;                              // sample size 
  int<lower = 2> N_g;                            // number of groups
  vector[N] y;                                   // response
  int<lower = 1, upper = N_g> group_id[N];       // group ID
}

transformed data{
  real y_mean;                                   // mean of y; see mu prior
  
  y_mean = mean(y); 
}

parameters {
  vector[2] mu;                                  // estimated group means and sd
  vector<lower = 0>[2] sigma;                    // Kruschke puts upper bound as well; ignored here
  real<lower = 0, upper = 100> nu;               // df for t distribution
}

model {
  // priors
  // note that there is a faster implementation of this for stan, 
  // and that the sd here is more informative than in Kruschke paper
  mu    ~ normal(y_mean, 10);                       
  sigma ~ cauchy(0, 5);
  
  // Based on Kruschke; makes average nu 29 
  // might consider upper bound, as if too large then might as well switch to normal
  nu    ~ exponential(1.0/29);                
  
  // likelihood
  for (n in 1:N) {
    y[n] ~ student_t(nu, mu[group_id[n]], sigma[group_id[n]]);
    
    // compare to normal; remove all nu specifications if you do this;
    //y[n] ~ normal(mu[group_id[n]], sigma[group_id[n]]);           
  }
}

generated quantities {
  vector[N] y_rep;                               // posterior predictive distribution
  real mu_diff;                                  // mean difference
  real cohens_d;                                 // effect size; see footnote 1 in Kruschke paper
  real CLES;                                     // common language effect size
  real CLES2;                                    // a more explicit approach; the mean should roughly equal CLES
  
  for (n in 1:N) {
    y_rep[n] = student_t_rng(nu, mu[group_id[n]], sigma[group_id[n]]);
  }
  
  mu_diff  = mu[1] - mu[2];
  cohens_d = mu_diff / sqrt(sum(sigma)/2);
  CLES     = normal_cdf(mu_diff / sqrt(sum(sigma)), 0, 1);
  CLES2    = student_t_rng(nu, mu[1], sigma[1]) - student_t_rng(nu, mu[2], sigma[2]) > 0;
}