test <- readMat("C:/Users/kaan/Documents/NatComm2023/MYELIN/DATA/MED.mat")
d_2<- data.frame( GS_std = standardize(test$MED[,1]) ,
                  ACW_std = standardize(test$MED[,2]) ,
                  MY_std = standardize(test$MED[,3]),
                  G = as.factor(test$MED[,5]))
# Synthetic data


mediated_1<- ulam( 
  alist( 
    GS_std ~ dnorm( mu , sigma ) , 
    mu <- a[G] + bACW[G] * ACW_std + bMY[G] * MY_std ,
    a[G] ~ dnorm( 0, 0.15 ) , 
    bACW[G] ~ dnorm( 0 , 0.27 ) ,
    bMY[G] ~ dnorm( 0 , 0.27 ) ,
    sigma ~ dexp( 1 ) ) , 
  data=d_2, chains = 4, cores = 4, log_lik = TRUE) 


mediated_2<- ulam( 
  alist( 
    GS_std ~ dnorm( mu , sigma ) , 
    # Interaction term is in the end
    mu <- a[G] + bACW[G] * ACW_std + bMY[G] * MY_std + bACWMY[G] * ACW_std * MY_std ,
    a[G] ~ dnorm( 0, 0.15 ) , 
    bACW[G] ~ dnorm( 0 , 0.25 ) ,
    bMY[G] ~ dnorm( 0 , 0.25 ) ,
    bACWMY[G] ~ dnorm( 0, 0.25 ),
    sigma ~ dexp( 1 ) ) , 
  data=d_2, chains = 4, cores = 4, log_lik = TRUE) 

######## 
# Basic model with GS is only related to myelin with multivariate priors and 
# partial pooling. Unfortunately, centered priors have divergent transition problems

# Synthetic data to check model performance and estimation
a <- 0 # Average gscorr content after standardized
b <- 0.5 #Average myelin difference after standardized
sigma_a <- 1 # std dev in intercepts 
sigma_b <- 0.5 # std dev in slopes 
rho <- (-0.7) # correlation between intercepts and slopes
Mu <- c( a , b )
# Covariance matrix
sigmas <- c(sigma_a,sigma_b) # standard deviations
Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix
# now matrix multiply to get covariance matrix 
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

n_regions <- 360

library(MASS) 
set.seed(5) # used to replicate example 
vary_effects <- mvrnorm( n_regions , Mu , Sigma )
a_G <- vary_effects[,1] 
bMY_G <- vary_effects[,2]
MY_std <- rnorm(n_regions, 0, 1)
G <- sample(c(rep(2, 33), rep(1, 327)))
mu <- a_G[G] + bMY_G[G] * MY_std 
sigma <- 0.5 
GS_std <- rnorm(n_regions, mu, sigma)

syn_data <- list(
  MY_std = MY_std,
  GS_std = GS_std,
  G = G
)

plot( a_G , bMY_G , col=rangi2 
      , xlab="intercepts (a_cafe)" , ylab="slopes (b_cafe)" )

# overlay population distribution 
library(ellipse) 
for ( l in c(0.1,0.3,0.5,0.8,0.99) ) 
  lines(ellipse(Sigma,centre=Mu,level=l),col=col.alpha("black",0.2))

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
    a ~ normal(0,0.5) , 
    bMY ~ normal(0,0.5) , 
    sigma_G ~ exponential(1), 
    sigma ~ exponential(1), 
    Rho ~ lkj_corr(2) ) , 
    data = syn_data, chains = 4, cores = 4, log_lik = TRUE) 

# Non - centered priors
m1_nc <- stan_model("C:/Users/kaan/Documents/NatComm2023/MYELIN/R/m1.stan")

fit_m1_nc <- sampling(
  m1_nc,
  data = syn_data ,
  thin = 4 ,
  cores = 4
)

# Prior and posterior correlation
post <- extract.samples( fit_m1_nc )
dens( post$Rho[,1,2] , xlim=c(-1,1) ) # posterior
R <- rlkjcorr( 1e4 , K=2 , eta=2 ) # prior
dens( R[,1,2] , add=TRUE , lty=2 )


# Posterior prediction plot
library(rethinking)
prior <- extract.samples( fit_m1_nc, n=1000 )

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

