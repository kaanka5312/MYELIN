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

d_subj <- list(MY_std = standardize(MY),
               GS_std = standardize(GS),
               subj = subj,
               G = G)


d_subj2 <- list(MY_std = standardize(MY) * (G-1),
                ACW_std = standardize(ACW0) * (G-1),
               GS_std = standardize(GS),
               G = subj)


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

N_cafes <- 10

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
N_visits <- 36
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


a <- 0 # average morning wait time 
b1 <- 0 # average difference afternoon wait time 
b2 <- 0

sigma_a <- 0.5 # std dev in intercepts 
sigma_b1 <- 0.25 # std dev in slopes 
sigma_b2 <- 0.25 # std dev in slopes 

#rho <- (-0.7) # correlation between intercepts and slopes
Mu <- c( a , b1, b2 )
# Covariance matrix
sigmas <- c( sigma_a,sigma_b1, sigma_b2 ) # standard deviations
Rho <- matrix(c(1, 0.7, 0.5, 
               0.7, 1, 0.35, 
               0.5, 0.35, 1), nrow = 3, ncol = 3)

# now matrix multiply to get covariance matrix 
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 10

library(MASS) 
set.seed(5) # used to replicate example 
vary_effects <- mvrnorm( N_cafes , Mu , Sigma )

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
    data = d_subj2, chains = 4, cores = 4, log_lik = TRUE, iter = 4000) 

# Non - centered priors
m1_nc <- stan_model("C:/Users/kaan/Documents/NatComm2023/MYELIN/R/m1.2.stan")

fit_m1_nc <- sampling(
  m1_nc,
  data = syn_data ,
  thin = 4 ,
  cores = 4
)

# ACW as a predictor added
m2 <-  ulam( 
  alist( 
    # Likelihood
    GS_std ~ dnorm( mu , sigma ) , 
    mu <- a_G[G] + bMY_G[G] * MY_std + bACW[G] * ACW_std,
    
    # Adaptive priors
    c(a_G,bMY_G)[G] ~ multi_normal( c(a,bMY), Rho, sigma_G) ,
    #vector[4]:c(a_G[G],bMY_G[G]) ~ multi_normal( c(a,bMY), Rho, sigma_G) ,
    #vector[2]:bMY_G[G] ~ multi_normal( 0, Rho, sigma_G) ,
    
    # Fixed priors
    a ~ normal(0,0.5) , 
    bMY ~ normal(-0.5,0.25) , 
    bACW[G] ~ dnorm(-0.5,0.25), # ACW added as a fixed prior
    sigma_ACW ~ exponential(1),
    sigma_G ~ exponential(1), 
    sigma ~ exponential(1), 
    Rho ~ lkj_corr(2)
    ) , 
  data = d_subj2, chains = 4, cores = 4, log_lik = TRUE, iter = 4000) 

# N O T E ! 
# There is nothing wrong with above multi - normal model. Only problem is 
# the following understanded from the syntehtic data: We are sampling only 2 
# samples from the multivariate dist. (As self and nonself), which hardens the 
# estimation of distrubition (because we only have 2 samples) ending up with 
# large SD posterior. Thus, another approach is taken. We apllied a varying 
# intercept and slope for each subject and take Beta as contrast between. As
# Can be understand from the synthetic data,non self regions are 0 and self's are 1.
# In below reflects varying intercepts for subject, however slopes are only dif-
# fered between self or nonself region groups. I think that is inferior.

# N O T E: Due to number of subjects, model is fitted in 45 mins.
m1.2 <-  ulam( 
  alist( 
    # Likelihood
    GS_std ~ dnorm( mu , sigma ) , 
    mu <- a[subj] + bMY[G] * MY_std ,
    
    # Adaptive priors
    a[subj] ~ dnorm( a_G, sigma_a ) ,
    bMY[G] ~ dnorm( bMY_G, sigma_b ) , 
    
    # Fixed priors
    a_G ~ normal(0,1) , 
    bMY_G ~ normal(0,1) , 
    sigma_a ~ exponential(1) , 
    sigma_b ~ exponential(1) , 
    sigma ~ exponential(1)
  ), data = d_subj, chains = 4, cores = 4, log_lik = TRUE, iter = 4000) 

# Non-centered version of m1.2 to avoid divergent errors

m1.2_nc <-  ulam( 
  alist( 
    # Likelihood
    GS_std ~ dnorm( mu , sigma ) , 
    #  Centered Version
    # mu <- a[subj] + bMY[G] * MY_std ,
    # Non - centered version
    mu <- a_G + z_a[subj] * sigma_a + bMY_G + z_b[G] * sigma_b, 
    
    # Adaptive priors
    # a[subj] ~ dnorm( a_G, sigma_a ) ,
    # bMY[G] ~ dnorm( bMY_G, sigma_b ) , 
    
    # Fixed priors
    a_G ~ normal(0,1) , 
    bMY_G ~ normal(0,1) ,
    
    z_a ~ normal(0,1),
    z_b ~ normal(0,1),
    
    sigma_a ~ exponential(1) , 
    sigma_b ~ exponential(1) , 
    sigma ~ exponential(1),
    
    # Generated Quantites
    gq> vector[subj]:a <- a_G + z_a * sigma_a,
    gq> vector[G]:bMY <- bMY_G + z_b * sigma_b  
      
  ), data = d_subj, chains = 4, cores = 4, log_lik = TRUE, iter = 4000) 


###### Plotting #######

# Posterior prediction plot
library(rethinking)
# Prior and posterior correlation
post <- extract.samples( m1 )

dens( post$Rho[,1,2] , xlim=c(-1,1) ) # posterior
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