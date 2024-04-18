a <- 3.5 # average morning wait time 
b <- (-1) # average difference afternoon wait time 
sigma_a <- 1 # std dev in intercepts 
sigma_b <- 0.5 # std dev in slopes 
rho <- (-0.7) # correlation between intercepts and slopes
Mu <- c( a , b )
# Covariance matrix
sigmas <- c(sigma_a,sigma_b) # standard deviations
Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix
# now matrix multiply to get covariance matrix 
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 20

library(MASS) 
set.seed(5) # used to replicate example 
vary_effects <- mvrnorm( N_cafes , Mu , Sigma )

a_cafe <- vary_effects[,1] 
b_cafe <- vary_effects[,2]

set.seed(22) 
N_visits <- 10 
afternoon <- rep(0:1,N_visits*N_cafes/2) 
#afternoon <- rnorm(200)
cafe_id <- rep( 1:N_cafes , each=N_visits ) 
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon 
sigma <- 0.5 # std dev within cafes 
wait <- rnorm( N_visits*N_cafes , mu , sigma ) 

d <- data.frame( cafe=cafe_id , afternoon=afternoon , wait=wait )

set.seed(867530) 
m14.1 <- ulam( alist( wait ~ normal( mu , sigma ), 
                      mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon,
                      c(a_cafe,b_cafe)[cafe] ~ multi_normal( c(a,b) , Rho , sigma_cafe ), 
                      a ~ normal(5,2), 
                      b ~ normal(-1,0.5), 
                      sigma_cafe ~ exponential(1), 
                      sigma ~ exponential(1), 
                      Rho ~ lkj_corr(2) ) , 
               data=d , chains=4 , cores=4,log_lik = TRUE )

m14.1nonc  <- ulam( alist( wait ~ normal( mu , sigma ), 
                                   mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon,
                                   #c(a_cafe,b_cafe)[cafe] ~ multi_normal( c(a,b) , Rho , sigma_cafe ), 
                                   vector[2]:c(a_cafe,b_cafe)[cafe] ~ multi_normal( c(a,b) , Rho , sigma_cafe ),
                                   a ~ normal(5,2), 
                                   b ~ normal(-1,0.5), 
                                   sigma_cafe ~ exponential(1), 
                                   sigma ~ exponential(1), 
                                   Rho ~ lkj_corr(2) ) , data=d , chains=4 , cores=4 )


post <- extract.samples(m14.1)
dens( post$Rho[,1,2] , xlim=c(-1,1) ) # posterior
R <- rlkjcorr( 1e4 , K=2 , eta=2 ) # prior
dens( R[,1,2] , add=TRUE , lty=2 )

# Testing in-home stan code
# Non - centered priors
stan_code<-'

data{
    vector[200] wait;
    int afternoon[200];
    int cafe[200];
}

parameters{
    real a;
    real b;
    vector<lower=0>[2] sigma_cafe;
    real<lower=0> sigma;
    cholesky_factor_corr[2] L_rho; // Means Rho matrix is a 2x2 matrix
    vector[2] z_alfa;
    vector[2] z_beta;
}

transformed parameters{
    vector[20] b_cafe;
    vector[20] a_cafe;
    a_cafe = (diag_pre_multiply(sigma_cafe, L_rho) * z_alfa);
    b_cafe = (diag_pre_multiply(sigma_cafe, L_rho) * z_beta);
  
}

model{
    vector[200] mu;
    L_rho ~ lkj_corr_cholesky( 2 );
    sigma ~ exponential( 1 );
    sigma_cafe ~ exponential( 1 );
    b ~ normal( -1 , 0.5 );
    a ~ normal( 5 , 2 );
    to_vector( z_alfa ) ~ normal( 0 , 1 );
    to_vector( z_beta ) ~ normal( 0 , 1 );
    
    for ( i in 1:200 ) {
        mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i];
    }
    wait ~ normal( mu , sigma );
}

generated quantities{
    vector[200] log_lik;
    vector[200] mu;
    matrix[2,2] Rho;
    Rho = multiply_lower_tri_self_transpose(L_rho);
    
 for ( i in 1:200 ) {
        mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i];
    }
    
    for ( i in 1:200 ) log_lik[i] = normal_lpdf( wait[i] | mu[i] , sigma );
}
'
stan_model_object <- stanc(model_code = stan_code)
model <- stan_model(stanc_ret = stan_model_object)
# Sample from the model
fit <- sampling(model, data = d, iter = 2000, chains = 4)

m14_nc <- stan_model("C:/Users/kaan/Documents/NatComm2023/MYELIN/R/m1_test.stan")

m14_nc_stan <- sampling(
  m14_nc,
  data = d ,
  thin = 4 ,
  cores = 4,
  chains = 4,
)
# compute unpooled estimates directly from data 
a1 <- sapply( 1:N_cafes , 
              function(i) mean(wait[cafe_id==i & afternoon==0]) ) 
b1 <- sapply( 1:N_cafes , 
              function(i) mean(wait[cafe_id==i & afternoon==1]) ) - a1

# extract posterior means of partially pooled estimates 
post <- extract.samples(m14.1) 
a2 <- apply( post$a_cafe , 2 , mean ) 
b2 <- apply( post$b_cafe , 2 , mean )

# plot both and connect with lines 
plot( a1 , b1 , xlab="intercept" , ylab="slope" , 
      pch=16 , col=rangi2 , ylim=c( min(b1)-0.1 , max(b1)+0.1 ) , 
      xlim=c( min(a1)-0.1 , max(a1)+0.1 ) ) 
points( a2 , b2 , pch=1 ) 
for ( i in 1:N_cafes ) lines( c(a1[i],a2[i]) , c(b1[i],b2[i]) )

# compute posterior mean bivariate Gaussian 
Mu_est <- c( mean(post$a) , mean(post$b) ) 
rho_est <- mean( post$Rho[,1,2] ) 
sa_est <- mean( post$sigma_cafe[,1] ) 
sb_est <- mean( post$sigma_cafe[,2] ) 
cov_ab <- sa_est*sb_est*rho_est 
Sigma_est <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )

######================ 14. 2 ====================#####
library(rethinking) 
data(chimpanzees) 
d <- chimpanzees 
d$block_id <- d$block 
d$treatment <- 1L + d$prosoc_left + 2L*d$condition

dat <- list( L = d$pulled_left, 
             tid = d$treatment, 
             actor = d$actor, 
             block_id = as.integer(d$block_id) )
set.seed(4387510) 
# Errors with divergent problems because we need non-centered priors
m14.2 <- ulam( alist( 
                      L ~ dbinom(1,p), 
                      logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid], 
                      # adaptive priors 
                      vector[4]:alpha[actor] ~ multi_normal(0,Rho_actor,sigma_actor), 
                      vector[4]:beta[block_id] ~ multi_normal(0,Rho_block,sigma_block), 
                      # fixed priors 
                      g[tid] ~ dnorm(0,1), 
                      sigma_actor ~ dexp(1), 
                      Rho_actor ~ dlkjcorr(4), 
                      sigma_block ~ dexp(1), 
                      Rho_block ~ dlkjcorr(4) ) , data=dat , chains=4 , cores=4 )

set.seed(4387510) 
m14.3 <- ulam( alist( L ~ binomial(1,p), 
                      logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
                      
                      # adaptive priors - non-centered
                      transpars> matrix[actor,4]:alpha <-
                        compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
                      transpars> matrix[block_id,4]:beta <-
                        compose_noncentered( sigma_block , L_Rho_block , z_block ),
                      matrix[4,actor]:z_actor ~ normal( 0 , 1 ),
                      matrix[4,block_id]:z_block ~ normal( 0 , 1 ),
                      
                      # fixed priors
                      g[tid] ~ normal(0,1),
                      vector[4]:sigma_actor ~ dexp(1),
                      cholesky_factor_corr[4]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
                      vector[4]:sigma_block ~ dexp(1),
                      cholesky_factor_corr[4]:L_Rho_block ~ lkj_corr_cholesky( 2 ),
                      
                      # compute ordinary correlation matrixes from Cholesky factors
                      gq> matrix[4,4]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
                      gq> matrix[4,4]:Rho_block <<- Chol_to_Corr(L_Rho_block)
) , data=dat , chains=4 , cores=4 , log_lik=TRUE )


m5.3_A <- quap( alist( 
  ## A -> D <- M 
  D ~ dnorm( mu , sigma ) ,
  mu <- a + bM*M + bA*A , 
  a ~ dnorm( 0 , 0.2 ) , 
  bM ~ dnorm( 0 , 0.5 ) , 
  bA ~ dnorm( 0 , 0.5 ) , 
  sigma ~ dexp( 1 ), 
  ## A -> M 
  M ~ dnorm( mu_M , sigma_M ), 
  mu_M <- aM + bAM*A, 
  aM ~ dnorm( 0 , 0.2 ), 
  bAM ~ dnorm( 0 , 0.5 ), 
  sigma_M ~ dexp( 1 ) ) , data = d )

post <- extract.samples( m5.3_A )
M_sim <- with( post , 
               sapply( 1:30 , 
                       function(i) rnorm( 1e3 , aM + bAM*A_seq[i] , sigma_M ) ) )


#######

a <- 3.5 # average morning wait time 
b <- (1) # average difference afternoon wait time 
c <- (1)

sigma_a <- 1 # std dev in intercepts 
sigma_b <- 0.5 # std dev in slopes 
sigma_c <- 0.5

#rho <- (0.7) # correlation between intercepts and slopes
Mu <- c( a , b, c)
# Covariance matrix

sigmas <- c(sigma_a,sigma_b, sigma_c) # standard deviations

Rho <- matrix(c(1, 0.7, 0.1, 
                0.7, 1, 0.07, 
                0.1, 0.07, 1), nrow = 3, ncol = 3)

#Rho <- matrix(c(1,rho,rho,1), ncol=2 )
# now matrix multiply to get covariance matrix 

Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 20

library(MASS) 
set.seed(5) # used to replicate example 
vary_effects <- mvrnorm( N_cafes , Mu , Sigma )

a_cafe <- vary_effects[,1] 
b_cafe <- vary_effects[,2]
c_cafe <- vary_effects[,3]

set.seed(22) 
N_visits <- 10
afternoon <- rep(0:1,N_visits*N_cafes/2) 
tired <- rnorm( N_visits*N_cafes, 1, 0.5  )

#afternoon <- rnorm(200)
cafe_id <- rep( 1:N_cafes , each=N_visits ) 
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon + c_cafe[cafe_id]*tired*afternoon
sigma <- 0.5 # std dev within cafes 
wait <- rnorm( N_visits*N_cafes , mu , sigma ) 

d <- data.frame( cafe=cafe_id , afternoon=afternoon , wait=wait )

set.seed(867530) 
m14.1 <- ulam( alist( wait ~ normal( mu , sigma ), 
                      mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon + c_cafe[cafe]*afternoon,
                      #c(a_cafe,b_cafe,c_cafe)[cafe] ~ multi_normal( c(a,b,c) , c(0.7,0.5,0.35) , sigma_cafe ), 
                      
                      # adaptive priors - non-centered
                      transpars> matrix[cafe,3]:coef <- compose_noncentered(, sigma_cafe, L_Rho, z_cafe ) ,
                      matrix[3,cafe]:z_cafe ~ normal( 0 , 1 ),
                      
                      gp > matrix[cafe,1]:a_cafe <- coef[cafe,1],
                      gp > matrix[cafe,1]:b_cafe <- coef[cafe,2],
                      gp > matrix[cafe,1]:c_cafe <- coef[cafe,3],

                      a ~ normal(5,2), 
                      b ~ normal(-1,0.5), 
                      c ~ normal(-1,0.5), 
                      vector[3]:sigma_cafe ~ exponential(1), 
                      sigma ~ exponential(1), 
                      cholesky_factor_corr[4]:L_Rho ~ ~ lkj_corr_cholesky( 2 ),
                      
                      gq> matrix[3,3]:Rho_cafe <- Chol_to_Corr(L_Rho)
                      ), 
               data=d , chains=4 , cores=4,log_lik = TRUE )

stan_code<-"
data{
    vector[2000] wait;
    int afternoon[2000];
    int cafe[2000];
}

parameters{
    real a;
    real b;
    real c;
    vector<lower=0>[3] sigma_cafe;
    real<lower=0> sigma;
    cholesky_factor_corr[3] L_rho; // Means Rho matrix is a 3x3 matrix
    vector[3] z[20];
}

transformed parameters{
    
    vector[20] c_cafe;
    vector[20] a_cafe;
    vector[20] b_cafe; 
    
    for ( j in 1:20 ){
    a_cafe[j] = a + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[1];
    b_cafe[j] = b + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[2];
    c_cafe[j] = c + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[3];
    }
}

model{
    vector[2000] mu;
    L_rho ~ lkj_corr_cholesky( 4 );
    sigma ~ exponential( 1 );
    sigma_cafe ~ exponential( 1 );
    c ~ normal( 1 , 0.5 );
    b ~ normal( 1 , 0.5 );
    a ~ normal( 5 , 2 );
    
    // Properly specify prior for each element of each vector in the z array
    for (i in 1:20) {
        z[i] ~ normal(0, 1);  // This applies the normal distribution to each 2D vector
    }
    
    for ( i in 1:2000 ) {
        mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i] + c_cafe[cafe[i]] * afternoon[i] ;
    }
    wait ~ normal( mu , sigma );
}

generated quantities{
    vector[2000] log_lik;
    vector[2000] mu;
    matrix[3,3] Rho;
    Rho = multiply_lower_tri_self_transpose(L_rho);
    
 for ( i in 1:2000 ) {
        mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i] + c_cafe[cafe[i]] * afternoon[i] ;
    }
    for ( i in 1:2000 ) log_lik[i] = normal_lpdf( wait[i] | mu[i] , sigma );
}

"




stan_model_object <- stanc(model_code = stan_code)
model <- stan_model(stanc_ret = stan_model_object)
# Sample from the model
fit <- sampling(model, data = d, iter = 2000, chains = 4)



library(MASS)

# Set parameters
set.seed(123)  # for reproducibility
N <- 100       # number of observations
alpha_true <- 0
beta1_true <- 1
beta2_true <- 2

# Stronger correlations for clearer effects
sigma <- 1
sigma_alpha <- 1
sigma_beta1 <- 1
sigma_beta2 <- 1
rho_alpha_beta1 <- 0.7
rho_alpha_beta2 <- 0.7
rho_beta1_beta2 <- 0.5
Covariance <- matrix(c(sigma_alpha^2, rho_alpha_beta1 * sigma_alpha * sigma_beta1, rho_alpha_beta2 * sigma_alpha * sigma_beta2,
                       rho_alpha_beta1 * sigma_alpha * sigma_beta1, sigma_beta1^2, rho_beta1_beta2 * sigma_beta1 * sigma_beta2,
                       rho_alpha_beta2 * sigma_alpha * sigma_beta2, rho_beta1_beta2 * sigma_beta1 * sigma_beta2, sigma_beta2^2),
                     nrow = 3, byrow = TRUE)

# Generate multivariate normal data for parameters
params <- mvrnorm(n = N, mu = c(alpha_true, beta1_true, beta2_true), Sigma = Covariance)

# Create covariates
X1 <- rnorm(N, mean = 2, sd = 1)
X2 <- rnorm(N, mean = 3, sd = 1.5)

# Generate outcome variable
mu <- params[, 1] + params[, 2] * X1 + params[, 3] * X2 
Y <- rnorm(N, mu, sigma)
# Prepare data for Stan
data_list <- list(N = N, Y = Y, X1 = X1, X2 = X2)




stan_gpt <- "
data {
  int<lower=0> N;
  vector[N] Y;
  vector[N] X1;
  vector[N] X2;
}

parameters {
  real a;
  real b;
  real c;
  vector[3] z[N];
  cholesky_factor_corr[3] L_rho;
  vector<lower=0>[3] sigma;
  real<lower=0> sigma_y;
}

transformed parameters {
  vector[N] beta_0;
  vector[N] beta_1;
  vector[N] beta_2;
  
  for ( j in 1:N ){
    beta_0[j] = a + (diag_pre_multiply(sigma, L_rho) * z[j])[1];
    beta_1[j] = b + (diag_pre_multiply(sigma, L_rho) * z[j])[2];
    beta_2[j] = c + (diag_pre_multiply(sigma, L_rho) * z[j])[3];
    }
}

model {
  L_rho ~ lkj_corr_cholesky(2);
  sigma ~ exponential( 1 );
  sigma_y ~ exponential( 1 );
  a ~ normal(0, 1);
  b ~ normal(1, 1);
  c ~ normal(2, 1);
  
  for (i in 1:N) {
        z[i] ~ normal(0, 1);  
    }
  
   // Correcting the model statement to use beta vectors
  for (i in 1:N) {
    Y[i] ~ normal(beta_0[i] + beta_1[i] * X1[i] + beta_2[i] * X2[i], sigma_y);
  }
}


generated quantities {
  matrix[3,3] Rho;
  Rho = multiply_lower_tri_self_transpose(L_rho);
}
"

stan_gpt <- "
data {
  int<lower=0> N;
  vector[N] Y;
  vector[N] X1;
  vector[N] X2;
}

parameters {
  real a;
  real b;
  real c;
  vector[3] z;
  cholesky_factor_corr[3] L_rho;
  vector<lower=0>[3] sigma;
  real<lower=0> sigma_y;
}

transformed parameters {
  vector[N] beta_0 = rep_vector(a, N) + rep_vector(sigma[1], N) .* (L_rho[1,1] * z[1] + L_rho[2,1] * z[2] + L_rho[3,1] * z[3]);
  vector[N] beta_1 = rep_vector(b, N) + rep_vector(sigma[2], N) .* (L_rho[1,2] * z[1] + L_rho[2,2] * z[2] + L_rho[3,2] * z[3]);
  vector[N] beta_2 = rep_vector(c, N) + rep_vector(sigma[3], N) .* (L_rho[1,3] * z[1] + L_rho[2,3] * z[2] + L_rho[3,3] * z[3]);
}

model {

  L_rho ~ lkj_corr_cholesky(2);
  sigma ~ exponential( 1 );
  sigma_y ~ exponential( 1 );
  a ~ normal(0, 1);
  b ~ normal(1, 1);
  c ~ normal(2, 1);
  z ~ normal(0, 1);  

  Y ~ normal(beta_0 + beta_1 .* X1 + beta_2 .* X2, sigma_y);
  
}


generated quantities {
  matrix[3,3] Rho;
  Rho = multiply_lower_tri_self_transpose(L_rho);
}
"


stan_model_object <- stanc(model_code = stan_gpt)
model <- stan_model(stanc_ret = stan_model_object)
# Sample from the model
fit.gpt <- sampling(model, data = data_list, iter = 2000, chains = 1, core  = 4)

m2 <-  ulam( 
  alist( 
    # Likelihood
    Y ~ dnorm( mu , sigma ) , 
    mu <- beta_0 + beta_1 * X1 + beta_2 * X2,
    beta_0 <- a + v[,1],
    beta_1 <- b + v[,2],
    beta_2 <- c + v[,3],
  
    # Adaptive priors
    vector[3]:v ~ multi_normal( 0, Rho, sigma_y) ,
    
    
    # Fixed priors
    a ~ normal(0,1) , 
    b ~ normal(1,1) , 
    c ~ dnorm(2,1), # ACW added as a fixed prior
    sigma ~ exponential(1),
    sigma_y ~ exponential(1), 
    cholesky_factor_corr[3]:L_Rho ~ lkj_corr_cholesky( 2 ),
    
    # Generated quantites
    gq> matrix[3,3]:Rho <<- Chol_to_Corr(L_Rho),
  ) , 
  data = data_list, chains = 4, cores = 4, log_lik = TRUE, iter = 4000) 
