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

xseq <- c(-2,2) 