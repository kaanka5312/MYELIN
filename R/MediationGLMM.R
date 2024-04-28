# Data preparation
# Averaged over subjects data
test <- readMat("C:/Users/kaan/Documents/NatComm2023/MYELIN/DATA/MED.mat")
d_2<- data.frame( GS_std = standardize(test$MED[,1]) ,
                  ACW_std = standardize(test$MED[,2]) ,
                  MY_std = standardize(test$MED[,3]),
                  G = as.factor(test$MED[,5]))

# All subjects and regions. Needed for varying intercept model
p1 <- readMat( "C:/Users/kaan/Documents/NatComm2023/MYELIN/DATA/INT_all.mat" )
p2 <- readMat( "C:/Users/kaan/Documents/NatComm2023/MYELIN/DATA/GSCORR.mat" )

ACW0 <- c(t(p1$ACW0.all))
MY <- c(t(p1$myelin.all))
GS <- c(p2$GSCORR.arr) # Careful that this in 360x100 formal already
subj <- rep(1:100, each = 360)
G <- rep(test$MED[,5],100)
self <- G - 1
nonself <- 1 - self

d_subj <- list(MY_std = standardize(MY),
               GS_std = standardize(GS),
               subj = subj,
               G = G)


d_subj2 <- list(MY_std = standardize(MY) * (G-1),
                ACW_std = standardize(ACW0) * (G-1),
               GS_std = standardize(GS),
               G = subj)


d_subj2 <- list(MY_std = standardize(MY) ,
                ACW_std = standardize(ACW0) ,
                GS_std = standardize(GS),
                G = subj,
                self = self ,
                nonself = nonself, 
                n_subj = 100,
                n_regions = 360
                )

######### Synthetic data ######
# Basic model with GS is only related to myelin with multivariate priors and 
# partial pooling. Unfortunately, centered priors have divergent transition problems
# Synthetic data to check model performance and estimation

# Synthetic Data 1 -- For only Myelin as causal
# Priors in the example of statistical rethinking chapter14
#a <- 3.5 # average morning wait time 
#b <- (-1) # average difference afternoon wait time 
#sigma_a <- 1 # std dev in intercepts 
#sigma_b <- 0.5 # std dev in slopes 

a <- 0 # average morning wait time 
b <- (-0.5) # average difference afternoon wait time 
sigma_a <- 0.5 # std dev in intercepts 
sigma_b <- 0.25 # std dev in slopes 


rho <- (-0.7) # correlation between intercepts and slopes
Mu <- c( a , b )
# Covariance matrix
sigmas <- c(sigma_a,sigma_b) # standard deviations
Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix
# now matrix multiply to get covariance matrix 
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 100

library(MASS) 
set.seed(5) # used to replicate example 
vary_effects <- mvrnorm( N_cafes , Mu , Sigma )

a_cafe <- vary_effects[,1] 
b_cafe <- vary_effects[,2]

set.seed(22) 
N_visits <- 360
#afternoon <- rep(0:1,N_visits*N_cafes/2)
afternoon <- rep(0:1,N_visits*N_cafes/2) * rnorm( N_visits*N_cafes )

cafe_id <- rep( 1:N_cafes , each=N_visits ) 
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon 
sigma <- 0.5 # std dev within cafes 
wait <- rnorm( N_visits*N_cafes , mu , sigma ) 

#d <- data.frame( cafe=cafe_id , afternoon=afternoon , wait=wait )
d_2 <- data.frame(G= cafe_id, MY_std = afternoon, GS_std = wait)


########## Synthetic data 2 ---- For Myelin + ACW
# ACW as a predictor added.
a <- 0 # average morning wait time 
b <- (-0.5) # average difference afternoon wait time 
sigma_a <- 0.5 # std dev in intercepts 
sigma_b <- 0.25 # std dev in slopes 


rho <- (-0.7) # correlation between intercepts and slopes
Mu <- c( a , b )
# Covariance matrix
sigmas <- c(sigma_a,sigma_b) # standard deviations
Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix
# now matrix multiply to get covariance matrix 
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 100

library(MASS) 
set.seed(5) # used to replicate example 
vary_effects <- mvrnorm( N_cafes , Mu , Sigma )

a_cafe <- vary_effects[,1] 
b_cafe <- vary_effects[,2]

set.seed(71)
b_ACW = c()
for (i in 1:N_cafes) {
  bACW_G = rnorm(1,mean = -0.5,sd = 0.25)
  z_ACW = rnorm(1)
  sigma_ACW = dexp(1)
  b_ACW[i] = bACW_G + z_ACW * sigma_ACW
}

set.seed(22) 
N_visits <- 360
#afternoon <- rep(0:1,N_visits*N_cafes/2)
#afternoon <- rep(0:1,N_visits*N_cafes/2) * rnorm( N_visits*N_cafes )

MY_std <- rep(0:1,N_visits*N_cafes/2) * rnorm( N_visits*N_cafes )
ACW_std <- rep(0:1,N_visits*N_cafes/2) * rnorm( N_visits*N_cafes )


cafe_id <- rep( 1:N_cafes , each=N_visits ) 
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*MY_std + b_ACW[cafe_id]*ACW_std
sigma <- 0.5 # std dev within cafes 
wait <- rnorm( N_visits*N_cafes , mu , sigma ) 

d_2 <- data.frame(G= cafe_id, MY_std = MY_std, ACW_std = ACW_std, GS_std = wait)

####

# Model captures the parameters in synthetic data well.
# Unfortunatelly centered priors are hard to converge
set.seed(4387510) 
m1 <-  ulam( 
  alist( 
    # Likelihood
    GS_std ~ dnorm( mu , sigma ) , 
    mu <- a_G[G] + bMY_G[G] * MY_std ,
    
    # Adaptive priors
    c(a_G,bMY_G)[G] ~ multi_normal( c(a,bMY), Rho, sigma_G) ,
    #vector[4]:c(a_G[G],bMY_G[G]) ~ multi_normal( c(a,bMY), Rho, sigma_G) ,
    #vector[2]:bMY_G[G] ~ multi_normal( 0, Rho, sigma_G) ,
    
    # Fixed priors
#    a ~ normal(5,2) , 
 #   bMY ~ normal(-1,0.5) , 
    a ~ normal(0,0.5) , 
    bMY ~ normal(-0.5,0.25) , 
    sigma_G ~ exponential(1), 
    sigma ~ exponential(1), 
    Rho ~ lkj_corr(2) ) , 
    data = d_2, chains = 4, cores = 4, log_lik = TRUE, iter = 4000) 

# Non - centered priors
m1.2 <- stan_model("C:/Users/kaan/Documents/NatComm2023/MYELIN/R/m1.2.stan")

fit_m1_nc <- sampling(
  m1.2,
  data = d_2 ,
  thin = 4 ,
  cores = 4,
  chains = 1
)

######### ACW as a predictor added ######
### Model example#####
## Synthetic Data ##
a <- 3.5   # Average morning wait time
b <- -1    # Average difference afternoon wait time
c <- 2     # Average effect of another covariate, like weekend

sigma_a <- 1   # Std dev in intercepts
sigma_b <- 0.5 # Std dev in slopes for b_cafe
sigma_c <- 0.5 # Std dev in slopes for c_cafe
rho <- -0.7    # Correlation between intercepts and slopes

Mu <- c(a, b, c)
sigmas <- c(sigma_a, sigma_b, sigma_c)
Rho <- matrix(c(1, rho, rho,
                rho, 1, 0,   # Assuming b_cafe and c_cafe are uncorrelated
                rho, 0, 1), nrow=3)  # Correlation matrix

# Now matrix multiply to get covariance matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 20
library(MASS)
set.seed(5)
vary_effects <- mvrnorm(N_cafes, Mu, Sigma)

a_cafe <- vary_effects[, 1]
b_cafe <- vary_effects[, 2]
c_cafe <- vary_effects[, 3]

set.seed(22)
N_visits <- 10
afternoon <- rep(0:1, N_visits * N_cafes / 2)
weekend <- sample(0:1, N_visits * N_cafes, replace = TRUE)  # New covariate, e.g., weekend or not

cafe_id <- rep(1:N_cafes, each = N_visits)
mu <- a_cafe[cafe_id] + b_cafe[cafe_id] * afternoon + c_cafe[cafe_id] * weekend
sigma <- 0.5
wait <- rnorm(N_visits * N_cafes, mu, sigma)

d <- data.frame(cafe = cafe_id, afternoon = afternoon, weekend = weekend, wait = wait)

### Stan Model ###
# Only MY -> GS <- ACW 
stan_code <-" 
data{
  vector[200] wait;
  int afternoon[200];
  int weekend[200];
  int cafe[200];
}

parameters{
  real a;
  real b;
  real c;
  vector<lower=0>[3] sigma_cafe;
  real<lower=0> sigma;
  cholesky_factor_corr[3] L_rho; // Means Rho matrix is a 2x2 matrix
  vector[3] z[20];
}

transformed parameters{
  vector[20] b_cafe;
  vector[20] a_cafe;
  vector[20] c_cafe;
  
  for ( j in 1:20 ){
    a_cafe[j] = a + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[1];
    b_cafe[j] = b + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[2];
    c_cafe[j] = c + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[3];
  }
}

model{
  vector[200] mu;
  L_rho ~ lkj_corr_cholesky( 2 );
  sigma ~ exponential( 1 );
  sigma_cafe ~ exponential( 1 );
  c ~ normal( 2 , 0.5 );
  b ~ normal( -1 , 0.5 );
  a ~ normal( 5 , 2 );
  
  // Properly specify prior for each element of each vector in the z array
  for (i in 1:20) {
    z[i] ~ normal(0, 1);  // This applies the normal distribution to each 2D vector
  }

  for ( i in 1:200 ) {
    mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i] + c_cafe[cafe[i]] * weekend[i];
  }
  wait ~ normal( mu , sigma );
}

generated quantities{
  vector[200] log_lik;
  vector[200] mu;
  matrix[3,3] Rho;
  Rho = multiply_lower_tri_self_transpose(L_rho);
  
  for ( i in 1:200 ) {
    mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i] + c_cafe[cafe[i]] * weekend[i];
  }
  for ( i in 1:200 ) log_lik[i] = normal_lpdf( wait[i] | mu[i] , sigma );
}
"
# This model has both:
# MY -> ACW
# MY -> GS <- ACW 
stan_code <-" 
data{
  vector[200] wait;
  int afternoon[200];
  int weekend[200];
  int cafe[200];
}

parameters{
  real a;
  real b;
  real c;
  vector<lower=0>[3] sigma_cafe;
  real<lower=0> sigma;
  cholesky_factor_corr[3] L_rho; // Means Rho matrix is a 2x2 matrix
  vector[3] z[20];
  real a_C;  // Intercept for modeling weekend as a function of afternoon
  real b_C;  // Slope for modeling weekend as a function of afternoon
}

transformed parameters{
  vector[20] b_cafe;
  vector[20] a_cafe;
  vector[20] c_cafe;
  
  for ( j in 1:20 ){
    a_cafe[j] = a + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[1];
    b_cafe[j] = b + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[2];
    c_cafe[j] = c + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[3];
  }
}

model{
  vector[200] mu;
  L_rho ~ lkj_corr_cholesky( 2 );
  sigma ~ exponential( 1 );
  sigma_cafe ~ exponential( 1 );
  c ~ normal( 2 , 0.5 );
  b ~ normal( -1 , 0.5 );
  a ~ normal( 5 , 2 );
  a_C ~ normal(0, 1);  // Reasonable prior for the intercept
  b_C ~ normal(0, 1);  // Reasonable prior for the slope
  
  // Properly specify prior for each element of each vector in the z array
  for (i in 1:20) {
    z[i] ~ normal(0, 1);  // This applies the normal distribution to each 2D vector
  }

 // Weekend is binary variable, thus bernoulli logit is used. This line consideres weekend
 // is related to afternoon and become mediator
  for ( i in 1: 200 ) {
   weekend[i] ~ bernoulli_logit(a_C + b_C * afternoon[i]);
  }

  for ( i in 1:200 ) {
    mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i] + c_cafe[cafe[i]] * weekend[i];
  }
  wait ~ normal( mu , sigma );
}

generated quantities{
  vector[200] log_lik;
  vector[200] mu;
  matrix[3,3] Rho;
  Rho = multiply_lower_tri_self_transpose(L_rho);
  
  for ( i in 1:200 ) {
    mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i] + c_cafe[cafe[i]] * weekend[i];
  }
  for ( i in 1:200 ) log_lik[i] = normal_lpdf( wait[i] | mu[i] , sigma );
}
"

stan_model_object <- stanc(model_code = stan_code)
model <- stan_model(stanc_ret = stan_model_object)
fit <- sampling(model, data = d, iter = 2000, chains = 4)

###########################################################
#################  T O Y  E X A M P L E  ##################
###########################################################

# This part creates several counterfactual plots in mediated model of cafes.
# Counterfactual 
posterior_samples <- extract(fit)
# Manual computation of the predicted outcomes
predict_wait_times <- function(samples, afternoon_scenario) {
  # Meaning across cafes 
  # Retrieve the number of posterior samples
  num_samples <- ncol(samples$a_cafe)
  N <- nrow(samples$a_cafe)
  
  # Compute the predictions for each sample
  predicted_wait <- matrix(ncol = num_samples, nrow = N)
  
    for (i in 1:num_samples) {
    # Calculate the mediator (weekend) under counterfactual afternoon scenario
    logit_p <- samples$a_C + samples$b_C * afternoon_scenario
    p_weekend <- 1 / (1 + exp(-logit_p))  # logistic function
    
    # Sample weekend occurrence
    weekend_scenario <- rbinom(N, 1, p_weekend)
    
    # Calculate the wait times
    predicted_wait[,i] <- samples$a_cafe[,i] + 
      samples$b_cafe[,i] * afternoon_scenario +
      samples$c_cafe[,i] * weekend_scenario
  }
  
  
  # Return the mean of predictions across all samples
  #apply(predicted_wait, 2, mean)
  return(predicted_wait)
}

# Generate counterfactual scenarios
afternoon_always <- rep(1, N)
afternoon_never <- rep(0, N)

# Generate predictions
cf_predictions_always <- predict_wait_times(posterior_samples, afternoon_always)
cf_predictions_never <- predict_wait_times(posterior_samples, afternoon_never)

#
mean_cafe_0 <- apply(cf_predictions_never, 2, mean)
mean_cafe_1 <- apply(cf_predictions_always, 2, mean)

# Shows difference between all afternoon or never, for each cafe
plot(mean_cafe_0)
points(mean_cafe_1,pch=16)
### Afternoon/ Myein set to 0. To simulate effect of weekend/ACW on W

# Parameters
n <- 200  # Number of elements in each array
num_arrays <- 1000  # Number of arrays to generate
p <- 0.5  # Probability of 1s in the arrays

# Generate the binary arrays
set.seed(123)  # Ensure reproducibility
binary_arrays <- replicate(num_arrays, sample(c(0, 1), n, replace = TRUE, prob = c(1 - p, p)))

# Calculate the sums of each array
array_sums <- colSums(binary_arrays)

predict_wait_times <- function(samples, afternoon_scenario, weekend_scenario) {
  # Meaning across cafes 
  # Retrieve the number of posterior samples
  num_samples <- ncol(samples$a_cafe)
  N <- nrow(samples$a_cafe)
  
  # Compute the predictions for each sample
  predicted_wait <- matrix(ncol = num_samples, nrow = N)
  
  for (i in 1:num_samples) {
    # Calculate the mediator (weekend) under counterfactual afternoon scenario
    logit_p <- samples$a_C + samples$b_C * afternoon_scenario
    p_weekend <- 1 / (1 + exp(-logit_p))  # logistic function
    
    # Sample weekend occurrence
    weekend_scenario <- rbinom(N, 1, p_weekend)
    
    # Calculate the wait times
    predicted_wait[,i] <- samples$a_cafe[,i] + 
      samples$b_cafe[,i] * afternoon_scenario +
      samples$c_cafe[,i] * weekend_scenario
  }
  
  
  # Return the mean of predictions across all samples
  apply(predicted_wait, 2, mean)
}

# Afternoon = 0. To delete the DAG from afternoon to both.
cf_predictions_never <- predict_wait_times(posterior_samples, afternoon_never)

col_list <- split(binary_arrays,col(binary_arrays))

mu_W <- sapply(1:1000, 
            function(i) predict_wait_times(posterior_samples, afternoon_never,
                                           binary_arrays[,i])  )
plot_mat <- rbind( colMeans(mu_W),
                   colSums(binary_arrays) )

ordered_plot_mat <- plot_mat[ ,order( plot_mat[2,] ) ]

plot(ordered_plot_mat[2,],
     ordered_plot_mat[1,])
abline(lm(ordered_plot_mat[1,]~
            ordered_plot_mat[2,]))
#############################################################
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=                

######## L I N E A R    D E C E N T E R E D ################
# Model that only makes MY -> GSCORR               
# Stan model ####
stan_linear <-"
data{
  vector[36000] GS_std;
  vector[36000] MY_std;
  int G[36000];
}
parameters{
  real a;
  real bMY;
  real bACW;
  vector<lower=0>[3] sigma_G;
  real<lower=0> sigma;
  cholesky_factor_corr[3] L_rho; // Means Rho matrix is a 3x3 matrix
  vector[3] z[100];
}

transformed parameters{
  vector[100] bMY_G;
  vector[100] a_G;
  for ( j in 1:100 ){
    
    a_G[j] = a + (diag_pre_multiply(sigma_G, L_rho) * z[j])[1];
    bMY_G[j] = bMY + (diag_pre_multiply(sigma_G, L_rho) * z[j])[2];
  } 
}
model{
  vector[36000] mu;
  L_rho ~ lkj_corr_cholesky( 2 );
  sigma ~ exponential( 1 );
  sigma_G ~ exponential( 1 );
  bMY ~ normal( -0.5 , 0.25 );
  a ~ normal( 0 , 0.5 );
  
  for (i in 1:100) {
    z[i] ~ normal(0, 1);  // This applies the normal distribution to each 2D vector
  }
  
  for ( i in 1:36000 ) {
    mu[i] = a_G[G[i]] + bMY_G[G[i]] * MY_std[i]  ;
  }
  GS_std ~ normal( mu , sigma );
}

generated quantities{
  vector[36000] log_lik;
  vector[36000] mu;
  matrix[3,3] Rho;
  Rho = multiply_lower_tri_self_transpose(L_rho);
  for ( i in 1:36000 ) {
    mu[i] = a_G[G[i]] + bMY_G[G[i]] * MY_std[i] ;
  }
  for ( i in 1:36000 ) log_lik[i] = normal_lpdf( GS_std[i] | mu[i] , sigma );
}
"
stan_model_object <- stanc(model_code = stan_linear)
model <- stan_model(stanc_ret = stan_model_object)
fit <- sampling(model, data = d, iter = 2000, chains = 4) #Synthetic
fit.linear <- sampling(model, data = d_subj2, iter = 2000, chains = 4) # Real


######## M E D I A T E D   D E C E N T E R E D ################               
# Model can capture the correlations in synthetic data 
#### Altered numbers and names, priors also #####
## Synthetic Data - Mediated ##
a <- 0 # average morning wait time 
b <- (-0.5) # average difference afternoon wait time 
c <- (0.5)

sigma_a <- 0.5   # Std dev in intercepts
sigma_b <- 0.25 # Std dev in slopes for b_cafe
sigma_c <- 0.25 # Std dev in slopes for c_cafe
rho <- -0.7    # Correlation between intercepts and slopes

Mu <- c(a, b, c)
sigmas <- c(sigma_a, sigma_b, sigma_c)
Rho_med <- matrix(c(1, rho, rho,
                rho, 1, 0,   # Assuming b_cafe and c_cafe are uncorrelated
                rho, 0, 1), nrow=3)  # Correlation matrix
Rho_dir <- matrix(c(1, rho, rho, 1), nrow=2)
# Now matrix multiply to get covariance matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 10
library(MASS)
set.seed(5)
vary_effects <- mvrnorm(N_cafes, Mu, Sigma)

a_cafe <- vary_effects[, 1]
b_cafe <- vary_effects[, 2]
c_cafe <- vary_effects[, 3]

set.seed(22)
N_visits <- 36
afternoon <- rep(0:1, N_visits * N_cafes / 2) * rnorm( N_visits*N_cafes )
# New covariate, e.g., weekend or not
weekend <- sample(0:1, N_visits * N_cafes, replace = TRUE) * rnorm( N_visits*N_cafes )  

cafe_id <- rep(1:N_cafes, each = N_visits)
mu <- a_cafe[cafe_id] + b_cafe[cafe_id] * afternoon + c_cafe[cafe_id] * weekend
sigma <- 0.5
wait <- rnorm(N_visits * N_cafes, mu, sigma)

d <- data.frame(G=cafe_id, MY_std = afternoon, ACW_std = weekend, GS_std = wait)
################## Synthetic data appropriate to the data ###################
#### N O T E  - I M P O R T A N T #####
# There is a problem with current version. We define data generation with 
# 0 and 1's with (G-1). However, model doesnt now that. This creates so variable
# Effect while MY -> ACW. Try with and without G-1. In order to solve this, 
# Next is done. Both self and nonself is defined with respeective variables 
# for each. This didn't delete, however shouldnt run See next attempt
## Synthetic Data - Mediated ##
a <- 0 # average GSCORR
ac <- 0 # average ACW between ACW-MY
b <- (-0.5) # average difference in myelin

bc <- 0.5 # Direct effect of MY_std on ACW_std
c <- (0.5) # average difference in ACW

sigma_a <- 0.5   # Std dev in intercepts
sigma_ac <- 0.5
sigma_b <- 0.25 # Std dev in slopes for b_cafe
sigma_c <- 0.25 # Std dev in slopes for c_cafe
sigma_bc <- 0.25

rho_med <- -0.7    # Correlation between intercepts and slopes
rho_dir <- 0.3

Mu_med <- c(a, b, c)
Mu_dir <- c(ac, bc)

sigmas_med <- c(sigma_a, sigma_b, sigma_c)
sigmas_dir <- c(sigma_ac, sigma_bc )

Rho_med <- matrix(c(1, rho_med, rho_med,
                    rho_med, 1, 0,   # Assuming b_cafe and c_cafe are uncorrelated
                    rho_med, 0, 1), nrow=3)  # Correlation matrix
Rho_dir <- matrix(c(1, rho_dir, rho_dir, 1), nrow=2)
# Now matrix multiply to get covariance matrix
Sigma_med <- diag(sigmas_med) %*% Rho_med %*% diag(sigmas_med)
Sigma_dir <- diag(sigmas_dir) %*% Rho_dir %*% diag(sigmas_dir)

n_subj <- 100
n_regions <- 360
G <- rep(test$MED[,5], n_subj)

library(MASS)
set.seed(5)
vary_effects <- mvrnorm(n_subj, Mu_med, Sigma_med)

a_subj <- vary_effects[, 1]
b_subj <- vary_effects[, 2]
c_subj <- vary_effects[, 3]

vary_effects <- mvrnorm(n_subj, Mu_dir, Sigma_dir)

ac_subj <- vary_effects[, 1]
bc_subj <- vary_effects[, 2]

subj_id <- rep(1:n_subj, each = n_regions)

# Direct
MY_std <- rnorm(n_subj*n_regions)
mu_ACW <- ac_subj[subj_id] + bc_subj[subj_id] * MY_std
sigma <- 0.5
ACW_std <- rnorm(n_subj * n_regions, mu_ACW, sigma)

# Includes direct effect of MY
mu <- a_subj[subj_id] + b_subj[subj_id] * MY_std + c_subj[subj_id] * ACW_std
# For MY -> ACW -> GS. Indirect
#mu <- a_subj[subj_id] + c_subj[subj_id] * ACW_std 
GS_std <- rnorm(n_subj * n_regions, mu, sigma)

d <- list(n_subj = n_subj, 
          n_regions = n_regions, 
          self = G-1,
          G=subj_id, 
          MY_std = MY_std * (G-1), 
          ACW_std = ACW_std * (G-1) , 
          GS_std = GS_std)

# Synthetic 
stan_syn <-"
data{
  int n_subj; // Total number of subjects
  int n_regions; //Total number of regions per subject
  vector[n_subj * n_regions ] self; // defines 1 as self 0 as nonself
  vector[n_subj * n_regions ] GS_std;
  vector[n_subj * n_regions] MY_std;
  vector[n_subj * n_regions] ACW_std;
  int G[n_subj * n_regions];
}
parameters{
  real a;
  real aC;
  real bMY;
  real bACW;
  real bC;
  vector<lower=0>[3] sigma_G;
  real<lower=0> sigma;
  vector<lower=0>[2] sigma_ACW;
  real<lower=0> sigma_2;
  cholesky_factor_corr[3] L_med; // Means Rho matrix is a 3x3 matrix
  cholesky_factor_corr[2] L_dir; // Means Rho matrix is a 2x2 matrix
  vector[3] z_med[n_subj];
  vector[2] z_dir[n_subj];
}

transformed parameters{
  vector[n_subj] bACW_G;
  vector[n_subj] bMY_G;
  vector[n_subj] a_G;
  for ( j in 1:n_subj ){
    a_G[j] = a + (diag_pre_multiply(sigma_G, L_med) * z_med[j])[1];
    bMY_G[j] = bMY + (diag_pre_multiply(sigma_G, L_med) * z_med[j])[2];
    bACW_G[j] = bACW + (diag_pre_multiply(sigma_G, L_med) * z_med[j])[3];
  } 
  
  vector[n_subj] bC_G;
  vector[n_subj] aC_G;
  for ( j in 1:n_subj ){
    aC_G[j] = aC + (diag_pre_multiply(sigma_ACW, L_dir) * z_dir[j])[1];
    bC_G[j] = bC + (diag_pre_multiply(sigma_ACW, L_dir) * z_dir[j])[2];
  } 
}

model{
  vector[n_subj * n_regions] mu;
  vector[n_subj * n_regions] mu_ACW;
  
  L_med ~ lkj_corr_cholesky( 2 );
  L_dir ~ lkj_corr_cholesky( 2 );
  
  sigma ~ exponential( 1 );
  sigma_G ~ exponential( 1 );
  sigma_ACW ~ exponential( 1 );
  sigma_2 ~ exponential( 1 );
  
  bACW ~ normal( 0.5 , 0.25 );
  bMY ~ normal( -0.5 , 0.25 );
  
  a ~ normal( 0 , 0.5 );
  aC ~ normal( -0.3, 0.25 ); // Direct effect of MY on ACW 
  bC ~ normal( -0.3 , 0.25 );
  
   
  for (i in 1:n_subj) {
    z_med[i] ~ normal(0, 1); 
    z_dir[i] ~ normal(0, 1);// This applies the normal distribution to each 2D vector
  }
  
    // MY -> ACW
    for ( i in 1:n_subj*n_regions ) {
     mu_ACW[i] = aC_G[G[i]] + bC_G[G[i]] * MY_std[i] * self[i];
    }
    
    ACW_std ~ normal( mu_ACW, sigma_2 );
    
    // MY -> GS <- ACW 
    for ( i in 1:n_subj*n_regions ) {
      mu[i] = a_G[G[i]] + bMY_G[G[i]] * MY_std[i] * self[i] + bACW_G[G[i]] * ACW_std[i] * self[i] ;
    }
    GS_std ~ normal( mu , sigma );
}

generated quantities{
  vector[n_subj*n_regions] mu;
  matrix[3,3] Rho_med;
  matrix[2,2] Rho_dir;
  Rho_med = multiply_lower_tri_self_transpose(L_med);
  Rho_dir = multiply_lower_tri_self_transpose(L_dir);
  for ( i in 1:n_subj*n_regions ) {
    mu[i] = a_G[G[i]] + bMY_G[G[i]] * MY_std[i] + bACW_G[G[i]] * ACW_std[i] ;
  }
}

"
stan_model_object <- stanc(model_code = stan_syn)
model <- stan_model(stanc_ret = stan_model_object)
#fit <- sampling(model, data = d, iter = 2000, chains = 4, cores= 4) #Synthetic
fit.mediated <- sampling(model, data = d_subj2, iter = 2000, chains = 4, cores= 6) # Real

post <- extract.samples(fit.mediated)
xseq <- seq(from = -2 , to = 2 , length = 100)
colors <- rainbow(10) # Generate distinct colors

# Counter-factual plot that showing relation between MY -> ACW, controlling MY               
plot(NULL,type = "l", ylim= c(-2.5,2.5),xlim=c(-2.5,2.5),
     xlab = "Standardized Myelin", ylab = "Standardized ACW", main= "Counterfactual Plot")
for (s in 1:10) {
  ACW_sim <- with(post, sapply(1:100, function(i) rnorm(1e3, aC_G[,s] + bC_G[,s]*xseq[i], sigma_2 )
  ))
  lines(xseq, y = colMeans(ACW_sim),col=col.alpha(colors[s],0.4)) 
  shade( apply(ACW_sim, 2, PI), xseq, col = col.alpha("gray",0.05))
}

# Total Effect of MY on GSCORR               
plot(NULL,type = "l", ylim= c(-3,3),xlim=c(-2.5,2.5),
     xlab = "Standardized Myelin", ylab = "Standardized GSCORR", main = "Total Effect of MY on GSCORR",sub= "Counterfactual Plot")
for (s in 1:10) {
ACW_sim <- with(post, sapply(1:100, function(i) rnorm(1e3, aC_G[,s] + bC_G[,s]*xseq[i], sigma_2 )
))
#lines(xseq, y = colMeans(ACW_sim)) 
#shade( apply(ACW_sim, 2, PI), xseq)

GS_sim <- sapply(1:100, function(i) rnorm(1e3, post$a_G[,s] + post$bMY_G[,s] * xseq[i] + 
                                            post$bACW_G[,s] * ACW_sim[,i],post$sigma))     

lines(xseq, y = colMeans(GS_sim), col=col.alpha(colors[s],0.4)) 
shade( apply(GS_sim, 2, PI), xseq,col = col.alpha("gray",0.05))
}

# Counter-factual Plot

plot(NULL,type = "l", ylim= c(-3,3),xlim=c(-2.5,2.5),
     xlab = "Standardized ACW", ylab = "Standardized GSCORR", main = "Counterfacted relationship")
for (s in 1:10) {
GS_sim <- sapply(1:100, function(i) rnorm(1e3, post$a_G[,s] + post$bMY_G[,s] * 0 + 
                                              post$bACW_G[,s] * xseq[i], post$sigma))     
lines(xseq, y = colMeans(GS_sim),col=col.alpha(colors[s],0.4)) 
shade( apply(GS_sim, 2, PI), xseq,col = col.alpha("gray",0.05))
}

#################### V A R Y I N G  E F F E C T   M E D I A T I O N ###########
#+-+--+--+--+--+--+- S E C O N D  A T T E M P T -+--+--+--+--+--+--+--+--+-
####### In direct pathway  #######

mean_parameters <- c( 0.2, -0.5, 0.2 , # Self mean 
                      -0.2, -0.3, 0.5 ) # Non-self mean

sigma_parameters <- c( 0.5, 0.25, 0.25,
                       0.5, 0.25, 0.25 )
names <- c("a1", "b1", "c1", "a2", "b2", "c2") #Row and colnames for convinience
Rho <- matrix(0, nrow = 6, ncol = 6)
colnames(Rho) <- names ; rownames(Rho) <- names
diag(Rho) <- 1

# Self parameters correlation 
Rho[1,2] <- (0.7) ; Rho[1,3] <- (0.7)
Rho[2,1] <- (0.7) ; Rho[3,1] <- (0.7)

# Non-self and self parameters correlation 
#Rho[1,4] <- (0.7) ; Rho[2,5] <- (0.7)
#Rho[4,1] <- (0.7) ; Rho[5,2] <- (0.7)

# Non self parameters intrinstic correlation
Rho[4,5] <- (0.5) ; Rho[4,6] <- (0.5)
Rho[5,4] <- (0.5) ; Rho[6,4] <- (0.5)

eigen(Rho)$values # Eigenvalues needs to be positive 
Sigma <- diag(sigma_parameters) %*% Rho %*% diag(sigma_parameters)

n_subj <- 10
n_regions <- 36

library(MASS)
set.seed(5)
vary_effects <- mvrnorm(n_subj, mean_parameters, Sigma)

# First column is self, second col non-self
a_subj_self <- vary_effects[, c(1)] 
b_subj_self <- vary_effects[, c(2)]
c_subj_self <- vary_effects[, c(3)]
a_subj_nonself <- vary_effects[, c(4)] 
b_subj_nonself <- vary_effects[, c(5)]
c_subj_nonself <- vary_effects[, c(6)]

######## Direct Pathway ########

mean_direct <- c(-0.3, 0.5, # Self
                 -0.1, 0.2 ) # Non-self

sigma_parameters <- c( 0.5, 0.25,
                       0.5, 0.25)

names <- c("a1", "b1", "a2", "b2") #Row and colnames for convinience
Rho <- matrix(0, nrow = 4, ncol = 4)
colnames(Rho) <- names ; rownames(Rho) <- names
diag(Rho) <- 1

# Self parameters correlation 
Rho[1,2] <- (0.3) ; Rho[1,3] <- (0.5)

# Non-self and self parameters correlation 
Rho[2,1] <- (0.3) ; Rho[3,1] <- (0.5)

# Non self parameters intrinstic correlation
Rho[4,3] <- (0.5) ; Rho[3,4] <- (0.3)


eigen(Rho)$values # Eigenvalues needs to be positive 
Sigma <- diag(sigma_parameters) %*% Rho %*% diag(sigma_parameters)

library(MASS)
set.seed(5)
vary_effects <- mvrnorm(n_subj, mean_direct, Sigma)

# First column is self, second col non-self
ac_subj_self <- vary_effects[, c(1)] 
bc_subj_self <- vary_effects[, c(2)]
ac_subj_nonself <- vary_effects[, c(3)] 
bc_subj_nonself <- vary_effects[, c(4)]


subj_id <- rep(1:n_subj, each = n_regions)
self <- rep(1:0, n_subj * n_regions/2)
nonself <- as.integer(!self)

# Direct
MY_std <- rnorm(n_subj*n_regions)

mu_ACW <- ac_subj_self[subj_id] * self +
  bc_subj_self[subj_id] * MY_std * self +
  ac_subj_nonself[subj_id] * nonself +
  bc_subj_nonself[subj_id] * MY_std * nonself

sigma <- 0.5
ACW_std <- rnorm(n_subj * n_regions, mu_ACW, sigma)

# Includes direct effect of MY
mu <- a_subj_self[subj_id] * self + 
  #      b_subj_self[subj_id] * MY_std * self+ 
  c_subj_self[subj_id] * ACW_std * self +
  a_subj_nonself[subj_id] * nonself + 
  #      b_subj_nonself[subj_id] * MY_std * nonself + 
  c_subj_nonself[subj_id] * ACW_std * nonself

# For MY -> ACW -> GS. Indirect
#mu <- a_subj[subj_id] + c_subj[subj_id] * ACW_std 
GS_std <- rnorm(n_subj * n_regions, mu, sigma)

d <- list(n_subj = n_subj, 
          n_regions = n_regions, 
          self = self,
          nonself = nonself,
          G=subj_id, 
          MY_std = MY_std , 
          ACW_std = ACW_std , 
          GS_std = GS_std)


# Synthetic 
stan_syn <-"
data{
  int n_subj; // Total number of subjects
  int n_regions; //Total number of regions per subject
  vector[n_subj * n_regions ] self; // defines 1 as self 0 as nonself
  vector[n_subj * n_regions ] nonself; // defines 1 as self 0 as nonself
  vector[n_subj * n_regions ] GS_std;
  vector[n_subj * n_regions] MY_std;
  vector[n_subj * n_regions] ACW_std;
  int G[n_subj * n_regions];
}
parameters{
  // MY -> GS <- ACW
  real a1;
  real b1;
  real c1;
  real a2;
  real b2;
  real c2;
  
  // MY -> ACW
  real ac1;
  real bc1;
  real ac2;
  real bc2;
  
  vector<lower=0>[6] sigma_G;
  real<lower=0> sigma;
  vector<lower=0>[4] sigma_ACW;
  real<lower=0> sigma_2;
  cholesky_factor_corr[6] L_med; // Means Rho matrix is a 6x6 matrix
  cholesky_factor_corr[4] L_dir; // Means Rho matrix is a 4x4 matrix
  vector[6] z_med[n_subj];
  vector[4] z_dir[n_subj];
}

transformed parameters{
  // Varying intercept self parameters
  vector[n_subj] a_subj_self;
  vector[n_subj] b_subj_self;
  vector[n_subj] c_subj_self;
  
  // Varying intercept nonself parameters
  vector[n_subj] a_subj_nonself;
  vector[n_subj] b_subj_nonself;
  vector[n_subj] c_subj_nonself;
  
  // Varying intercept for MY -> ACW
   vector[n_subj] ac_subj_self;
   vector[n_subj] bc_subj_self;
   vector[n_subj] ac_subj_nonself;
   vector[n_subj] bc_subj_nonself;
  
  // MY -> GS <- ACW
  for ( j in 1:n_subj ){
    a_subj_self[j] = a1 + (diag_pre_multiply(sigma_G, L_med) * z_med[j])[1];
    b_subj_self[j] = b1 + (diag_pre_multiply(sigma_G, L_med) * z_med[j])[2];
    c_subj_self[j] = c1 + (diag_pre_multiply(sigma_G, L_med) * z_med[j])[3];
    a_subj_nonself[j] = a2 + (diag_pre_multiply(sigma_G, L_med) * z_med[j])[4];
    b_subj_nonself[j] = b2 + (diag_pre_multiply(sigma_G, L_med) * z_med[j])[5];
    c_subj_nonself[j] = c2 + (diag_pre_multiply(sigma_G, L_med) * z_med[j])[6];
  } 
  
  // MY -> ACW
  for ( j in 1:n_subj ){
    ac_subj_self[j] = ac1 + (diag_pre_multiply(sigma_ACW, L_dir) * z_dir[j])[1];
    bc_subj_self[j] = bc1 + (diag_pre_multiply(sigma_ACW, L_dir) * z_dir[j])[2];
    ac_subj_nonself[j] = ac2 + (diag_pre_multiply(sigma_ACW, L_dir) * z_dir[j])[3];
    bc_subj_nonself[j] = bc2 + (diag_pre_multiply(sigma_ACW, L_dir) * z_dir[j])[4];
  } 
}

model{
  // Hyperpriors
  a1 ~ normal( 0, 0.5);
  b1 ~ normal( 0, 0.5);
  c1 ~ normal( 0, 0.5);
  a2 ~ normal( 0, 0.5);
  b2 ~ normal( 0, 0.5);
  c2 ~ normal( 0, 0.5);
  
  ac1 ~ normal( 0, 0.5 );
  bc1 ~ normal( 0, 0.5 );
  ac2 ~ normal( 0, 0.5 );
  bc2 ~ normal( 0, 0.5 );
  
  vector[n_subj * n_regions] mu;
  vector[n_subj * n_regions] mu_ACW;
  
  L_med ~ lkj_corr_cholesky( 2 );
  L_dir ~ lkj_corr_cholesky( 2 );
  
  sigma ~ exponential( 1 );
  sigma_G ~ exponential( 1 );
  sigma_ACW ~ exponential( 1 );
  sigma_2 ~ exponential( 1 );
  
  for (i in 1:n_subj) {
    z_med[i] ~ normal(0, 1); 
    z_dir[i] ~ normal(0, 1);// This applies the normal distribution to each 2D vector
  }
  
    // MY -> ACW
    for ( i in 1:n_subj*n_regions ) {
     mu_ACW[i] = ac_subj_self[G[i]] * self[i] + 
                 bc_subj_self[G[i]] * MY_std[i] * self[i] +
                 ac_subj_nonself[G[i]] * nonself[i] + 
                 bc_subj_nonself[G[i]] * MY_std[i] * nonself[i] ;
    }
    
    ACW_std ~ normal( mu_ACW, sigma_2 );
    
    // MY -> GS <- ACW 
    for ( i in 1:n_subj*n_regions ) {
      mu[i] = a_subj_self[G[i]] * self[i] +
              b_subj_self[G[i]] * MY_std[i] * self[i] + 
              c_subj_self[G[i]] * ACW_std[i] * self[i] +
              a_subj_nonself[G[i]] * nonself[i] +
              b_subj_nonself[G[i]] * MY_std[i] * nonself[i] + 
              c_subj_nonself[G[i]] * ACW_std[i] * nonself[i] ;
    }
    GS_std ~ normal( mu , sigma );
}

generated quantities{
  matrix[6,6] Rho_med;
  matrix[4,4] Rho_dir;
  Rho_med = multiply_lower_tri_self_transpose(L_med);
  Rho_dir = multiply_lower_tri_self_transpose(L_dir);
  
}

"
stan_model_object <- stanc(model_code = stan_syn)
model <- stan_model(stanc_ret = stan_model_object)
#fit <- sampling(model, data = d, iter = 2000, chains = 4, cores= 4) #Synthetic
fit.mediated <- sampling(model, data = d_subj2, iter = 2000, chains = 4, cores= 6) # Real
save(fit.mediated, file = "/media/kaansocat/Elements/EIB/DATA/MultMed.RData" )

" 
  vector[n_subj * n_regions] GS_std_prior;
  vector[n_subj * n_regions] ACW_std_prior;
// prior predictive check
  for ( i in 1:n_subj*n_regions ) {
    // Generate ACW_std based on prior distributions
    ACW_std_prior[i] = normal_rng(
      ac1 * self[i] + bc1 * MY_std[i] * self[i] + 
      ac2 * nonself[i] + bc2 * MY_std[i] * nonself[i],
      sigma_2
    );
  
   // Generate GS_std based on prior distributions and using ACW_std_prior
    GS_std_prior[i] = normal_rng(
      a1 * self[i] + b1 * MY_std[i] * self[i] + c1 * ACW_std_prior[i] * self[i] + 
      a2 * nonself[i] + b2 * MY_std[i] * nonself[i] + c2 * ACW_std_prior[i] * nonself[i],
      sigma
    );
  } "
post <- extract.samples(fit)
xseq <- seq(from = -2 , to = 2 , length = 100)
colors <- rainbow(2) # Generate distinct colors for self and nonself

####### Counter-factual plot that showing relation between MY -> ACW, controlling MY #######               
# Meaning across subjects for visual purposes
n_subj = 100
colMeanPlot <- matrix(NA, ncol = n_subj, nrow = length(xseq))

plot(NULL,type = "l", ylim= c(-2.5,2.5),xlim=c(-2.5,2.5),
     xlab = "Standardized Myelin", ylab = "Standardized ACW", main= "Counterfactual Plot")

for (s in 1:n_subj) {
  ACW_sim <- with(post, sapply(1:100, function(i) rnorm(1e3, 
                                                        ac_subj_nonself[,s] + bc_subj_nonself[,s]*xseq[i], 
                                                        sigma_2 )
  ))
  #  lines(xseq, y = colMeans(ACW_sim),col=col.alpha(colors[1],0.4)) 
  #  shade( apply(ACW_sim, 2, PI), xseq, col = col.alpha(colors[1],0.05))
  colMeanPlot[,s] = colMeans(ACW_sim)
  
}
lines(xseq, y = rowMeans(colMeanPlot), lwd=2, col = col.alpha(colors[1],0.8))

colMeanPlot <- matrix(NA, ncol = n_subj, nrow = length(xseq))
for (s in 1:n_subj) {
  ACW_sim <- with(post, sapply(1:100, function(i) rnorm(1e3, 
                                                        ac_subj_self[,s] + bc_subj_self[,s]*xseq[i], 
                                                        sigma_2 )
  ))
  #  lines(xseq, y = colMeans(ACW_sim),col=col.alpha(colors[2],0.4)) 
  #  shade( apply(ACW_sim, 2, PI), xseq, col = col.alpha(colors[2],0.05))
  colMeanPlot[,s] = colMeans(ACW_sim)
}
lines(xseq, y = rowMeans(colMeanPlot), lwd=2, col = col.alpha(colors[2],0.8))


############### Total Effect of MY on GSCORR ###################

#=+=+=+=+=+=+=+=+ N O N - S E L F =+=+=+=+=+=+=+=+=+=+ #

colMeanPlot <- matrix(NA, ncol = n_subj, nrow = length(xseq))

plot(NULL,type = "l", ylim= c(-3,3),xlim=c(-2.5,2.5),
     xlab = "Standardized Myelin", ylab = "Standardized GSCORR", main = "Total Effect of MY on GSCORR",sub= "Counterfactual Plot")

for (s in 1:n_subj) {
  ACW_sim <- with(post, sapply(1:100, function(i) rnorm(1e3, 
                                                        ac_subj_nonself[,s] + bc_subj_nonself[,s]*xseq[i],  
                                                        sigma_2 )
  ))
  #lines(xseq, y = colMeans(ACW_sim)) 
  #shade( apply(ACW_sim, 2, PI), xseq)
  
  GS_sim <- with(post, sapply(1:100, function(i) rnorm(1e3, 
                                                       a_subj_nonself[,s] + b_subj_nonself[,s] * xseq[i] + 
                                                         c_subj_nonself[,s] * ACW_sim[,i],
                                                       sigma)) )     
  colMeanPlot[,s] = colMeans(GS_sim)
  #lines(xseq, y = colMeans(GS_sim), col=col.alpha(colors[1],0.4)) 
  #shade( apply(GS_sim, 2, PI), xseq,col = col.alpha(colors[1],0.05))
}
lines(xseq, y = rowMeans(colMeanPlot), lwd=2, col = col.alpha(colors[1],0.8))

#=+=+=+=+=+=+=+=+ S E L F =+=+=+=+=+=+=+=+=+=+ #

colMeanPlot <- matrix(NA, ncol = n_subj, nrow = length(xseq))
for (s in 1:n_subj) {
  ACW_sim <- with(post, sapply(1:100, function(i) rnorm(1e3, 
                                                        ac_subj_self[,s] + bc_subj_self[,s]*xseq[i],  
                                                        sigma_2 )
  ))
  #lines(xseq, y = colMeans(ACW_sim)) 
  #shade( apply(ACW_sim, 2, PI), xseq)
  
  GS_sim <- with(post, sapply(1:100, function(i) rnorm(1e3, 
                                                       a_subj_self[,s] + b_subj_self[,s] * xseq[i] + 
                                                         c_subj_self[,s] * ACW_sim[,i],
                                                       sigma)) )     
  colMeanPlot[,s] = colMeans(GS_sim)
  #lines(xseq, y = colMeans(GS_sim), col=col.alpha(colors[2],0.4)) 
  #shade( apply(GS_sim, 2, PI), xseq,col = col.alpha(colors[2],0.05))
}
lines(xseq, y = rowMeans(colMeanPlot), lwd=2, col = col.alpha(colors[2],0.8))

############################## Counter-factual Plot ########################

#=+=+=+=+=+=+=+=+ N O N - S E L F =+=+=+=+=+=+=+=+=+=+ #

colMeanPlot <- matrix(NA, ncol = n_subj, nrow = length(xseq))

plot(NULL,type = "l", ylim= c(-3,3),xlim=c(-2.5,2.5),
     xlab = "Standardized ACW", ylab = "Standardized GSCORR", main = "Counterfacted relationship")
for (s in 1:n_subj) {
  GS_sim <- with(post,sapply(1:100, function(i) rnorm(1e3, a_subj_nonself[,s] + 
                                                        b_subj_nonself[,s] * 0 + 
                                                        c_subj_nonself[,s] * xseq[i], 
                                                      sigma)) )    
  colMeanPlot[,s] = colMeans(GS_sim)
  # lines(xseq, y = colMeans(GS_sim),col=col.alpha(colors[1],0.4)) 
  #  shade( apply(GS_sim, 2, PI), xseq,col = col.alpha(colors[1],0.05))
}
lines(xseq, y = rowMeans(colMeanPlot), lwd=2, col = col.alpha(colors[1],0.8))

#=+=+=+=+=+=+=+=+ S E L F =+=+=+=+=+=+=+=+=+=+ #

colMeanPlot <- matrix(NA, ncol = n_subj, nrow = length(xseq))
for (s in 1:n_subj) {
  GS_sim <- with(post,sapply(1:100, function(i) rnorm(1e3, a_subj_self[,s] + 
                                                        b_subj_self[,s] * 0 + 
                                                        c_subj_self[,s] * xseq[i], 
                                                      sigma)) )    
  colMeanPlot[,s] = colMeans(GS_sim)
  #lines(xseq, y = colMeans(GS_sim),col=col.alpha(colors[2],0.4)) 
  #shade( apply(GS_sim, 2, PI), xseq,col = col.alpha(colors[2],0.05))
}
lines(xseq, y = rowMeans(colMeanPlot), lwd=2, col = col.alpha(colors[2],0.8))

############## PRIOR PREDICTIVE CHECK ################
a1 <- rnorm(1e3, 0, 0.5 )
b1 <- rnorm(1e3, 0, 0.5 )
c1 <- rnorm(1e3, 0, 0.5 )

a2 <- rnorm(1e3, 0, 0.5 )
b2 <- rnorm(1e3, 0, 0.5 )
c2 <- rnorm(1e3, 0, 0.5)

ac1 <- rnorm(1e3, 0, 0.5)
bc1 <- rnorm(1e3, 0, 0.5)
ac2 <- rnorm(1e3, 0, 0.5)
bc2 <- rnorm(1e3, 0, 0.5)

self <- rbinom(1000, 1, 1)  # Binary, assuming about half are 'self'
nonself <- 1 - self  # Directly opposite of self
MY_std <- rnorm(1000)  # Standard normal distribution for MY_std
sigma_2 <- 0.5  # Assuming a constant standard deviation for simplicity
sigma <- 0.5

xseq <- c(-2,2)
# Compute the model for ACW_std_prior based on the parameters and data
# Simulate ACW_std_prior from the normal distribution with mu_ACW and sigma_2
ACW_sim <- sapply(1:2, function(i) rnorm(1e3,ac1 * self + bc1 * xseq[i] * self +
                                           ac2 * nonself + bc2 * xseq[i] * nonself, sigma_2 ) )
plot(NULL,type = "l", ylim= c(-3,3),xlim=c(-2.5,2.5),
     xlab = "Standardized ACW", ylab = "Standardized GSCORR", main = "Counterfacted relationship")
#lines(xseq, y = colMeans(ACW_sim),col=col.alpha(colors[1],1), lwd=2) 
for (i in 1:200) lines(xseq, ACW_sim[i,],col=col.alpha("black",0.4), lwd=2 )

# Calculate mu_GS based on the model specified

GS_sim <- sapply(1:2, function(i) rnorm(1e3,mean = a1 * self + b1 * xseq[i] * self + c1 * ACW_sim[,i] * self +
                                          a2 * nonself + b2 * xseq[i] * nonself + c2 * ACW_sim[,i] * nonself,
                                        sigma) )      

plot(NULL,type = "l", ylim= c(-3,3),xlim=c(-2.5,2.5),
     xlab = "Standardized ACW", ylab = "Standardized GSCORR", main = "Counterfacted relationship")
#lines(xseq, y = colMeans(ACW_sim),col=col.alpha(colors[1],1), lwd=2) 
for (i in 1:200) lines(xseq, GS_sim[i,],col=col.alpha("black",0.4), lwd=2 )

# Define the number of simulations and subjects
n <- 1000
n_subj <- 10  # Assuming 10 subjects for simplicity

# Generate parameters using normal distributions based on priors
a_subj_self <- rnorm(n_subj, 0, 0.5)
b_subj_self <- rnorm(n_subj, -0.5, 0.25)
c_subj_self <- rnorm(n_subj, 0, 0.25)

# Define the standard deviation for the outcomes
sigma <- 0.25

# Create a sequence over which to evaluate the model
xseq <- seq(-1, 1, length.out = 100)
# Pre-allocate a matrix to store simulation results
GS_sim <- matrix(nrow = n, ncol = length(xseq))

# Simulate outcomes for each value in xseq
for (i in seq_along(xseq)) {
  GS_sim[, i] <- rnorm(n, 
                       mean = a_subj_self + c_subj_self * xseq[i],
                       sd = sigma)
}

plot(NULL,type = "l", ylim= c(-2,2),xlim=c(-2,2),
     xlab = "Standardized ACW", ylab = "Standardized GSCORR", main = "Counterfacted relationship")
lines(xseq, y = colMeans(GS_sim),col=col.alpha(colors[2],1)) 
shade( apply(GS_sim, 2, PI), xseq,col = col.alpha(colors[2],0.2))
###### Plotting #######

# Posterior prediction plot
library(rethinking)
# Prior and posterior correlation
post <- extract.samples( m1 )

dens( post$Rho_med[,1,2] , xlim=c(-1,1) ) # posterior
R <- rlkjcorr( 1e4 , K=2 , eta=2 ) # prior
dens( R[,1,2] , add=TRUE , lty=2 )

dens( post$bMY_G[,1] , xlim=c(-1,1) ) # posterior
prior<- mvrnorm( 1e4 , Mu , Sigma ) 
dens( prior[,1] , add=TRUE , lty=2 )
prior <- extract.samples( m1, n=1000 )

N <- 1000 # total number of observations
N_g <- 2 # number of groups
group_id <- sample(1:N_g, N, replace = TRUE) # Randomly assign observations to 2 groups for demonstration
xseq <- seq(-2,2, length.out=30) 

y_simulated <- matrix(nrow = N, ncol = length(xseq)+1) # First column is for group

y_sim1 <- with(prior,
            sapply(1:30, 
            function(i) rnorm(1e3, a_G[1] + bMY_G[1] * xseq[i], sigma)))

y_sim2 <- with(prior,
               sapply(1:30, 
               function(i) rnorm(1e3, a_G[2] + bMY_G[2] * xseq[i], sigma)))

sim_list = list(y_sim1,y_sim2)

# Wrong plotting bu just here for if structure needed
# for (i in 1:N) {
#   group_id_n <- sample(1:N_g, 1) # Randomly assign each data point to a group
#   # and saving to the first column
#   for (n in 1:length(xseq)) {
#     
#     y_simulated[i,1] <- group_id_n # Saving the which group belonged to
#     
#     # Parameters
#     a_Gi <- prior$a_G[i,group_id_n]
#     bMY_Gi <- prior$bMY_G[i,group_id_n]
#     mu_i <- a_Gi + bMY_G[G] * xseq[n]
#     sigma_i <- prior$sigma[i]
#     y_simulated[i, n+1] <- rnorm(1, mean = mu_i, sd= sigma_i) # Generate the data point
#     
#   }
# }


# Plotting 
rows_to_plot <- seq(2, 1000, by = 20) # For example, plotting every 200th row
colors <- rainbow(N_g) # Generate distinct colors
alpha = 0.2

# Initialize the plot with the first density
dens <- density(y_simulated[rows_to_plot[1], ])
plot(dens, main = "Density of Simulated Data", xlab = "Simulated Values", ylab = "Density", type = 'l', col = 1)

# Add the rest of the densities using lines()
for (i in 1:length(rows_to_plot)) {
  dens <- density(y_simulated[rows_to_plot[i], ])
  lines(dens, col = col.alpha(colors[y_simulated[i,1]],0.4)) # Use different colors 
}

dens <- density(sim_list[[1]][1, ])
plot(dens, main = "Density of Simulated Data", xlab = "Simulated Values", ylab = "Density", type = 'l', col = 1)

# Add the rest of the densities using lines()
for (i in 1:length(xseq)) {
  dens <- density(sim_list[[1]][i, ])
  lines(dens, col = col.alpha(colors[1],0.4)) # Use different colors 
  
  dens <- density(sim_list[[2]][i, ])
  lines(dens, col = col.alpha(colors[2],0.4)) # Use different colors 
  
}

# Mean lines and shade
mu_mean1 <- apply(y_sim1,2,mean) 
mu_PI1 <- apply(y_sim1,2,PI) 
mu_mean2 <- apply(y_sim2,2,mean) 
mu_PI2 <- apply(y_sim2,2,PI) 

plot( GS_std ~ MY_std  , data=d_2[d_2$G==1,],
      xlim = c(-2,2), ylim = c(-2,2), col = col.alpha(colors[1],0.3) ) 
lines( xseq , mu_mean1 , lwd=2, col=colors[1] ) 
shade( mu_PI1 , xseq , col = col.alpha(colors[1],0.15))
points( GS_std ~ MY_std  , data=d_2[d_2$G==2,], col = col.alpha(colors[2],0.3)  ) 
lines( xseq , mu_mean2 , lwd=2, col=colors[2] ) 
shade( mu_PI2 , xseq , col = col.alpha(colors[2],0.15))

## New 

# compute posterior mean bivariate Gaussian 
post <- extract.samples(m14.1)

Mu_est <- c( mean(post$a) , mean(post$b) ) 
rho_est <- mean( post$Rho[,1,2] ) 
sa_est <- mean( post$sigma_cafe[,1] ) 
sb_est <- mean( post$sigma_cafe[,2] ) 
cov_ab <- sa_est*sb_est*rho_est 
Sigma_est <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )

post <- extract.samples(m1)

a2 <- apply( post$a_G , 2 , mean ) 
b2 <- apply( post$bMY_G , 2 , mean )

Mu_est <- c( mean(post$a) , mean(post$bMY) ) 
rho_est <- mean( post$Rho[,1,2] ) 
sa_est <- mean( post$sigma_G[,1] ) 
sb_est <- mean( post$sigma_G[,2] ) 
cov_ab <- sa_est*sb_est*rho_est 
Sigma_est <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )

plot(a2,b2, xlab = "intercept", ylab = "slope", pch = 16, 
     col=rangi2 , ylim=c( min(b2)-0.1 , max(b2)+0.1 ) , xlim=c( min(a2)-0.1 , max(a2)+0.1) )
library(ellipse)
for ( l in c(0.1,0.3,0.5,0.8,0.99) ) lines(ellipse(Sigma_est,centre=Mu_est,level=l), col=col.alpha("black",0.2),lwd=2)

# Convert varying effects to GS
GS_ns <- a2
GS_s <- a2 + b2

plot( GS_ns, GS_s, xlab = "standardized GS nonself", ylab = "standardized GS Self",
      pch = 16, col = rangi2, 
      ylim = c( min(GS_s)-0.1 , max(GS_s)+0.1 ),
      xlim=c( min(GS_ns)-0.1 , max(GS_ns)+0.1 ) )
abline( a=0 , b=1 , lty=2, lwd=2 )

# now shrinkage distribution by simulation 
v <- mvrnorm( 1e4 , Mu_est , Sigma_est ) 
v[,2] <- v[,1] + v[,2] # calculate afternoon wait 
Sigma_est2 <- cov(v) 
Mu_est2 <- Mu_est 
Mu_est2[2] <- Mu_est[1]+Mu_est[2]

library(ellipse) 
for ( l in c(0.1,0.3,0.5,0.8,0.99) ) lines(ellipse(Sigma_est2,centre=Mu_est2,level=l), col=col.alpha("black",0.5), lwd=2)

