data{
    vector[360] GS_std;
    vector[360] MY_std;
    int G[360];
}
parameters{
    real a;
    real bMY;
    vector<lower=0>[2] sigma_G;
    real<lower=0> sigma;
    cholesky_factor_corr[2] L_rho; // Means Rho matrix is a 2x2 matrix
    vector[2] z[2];
}

transformed parameters{
    vector[2] bMY_G;
    vector[2] a_G;
    for ( j in 1:2 ){
       
    a_G[j] = a + (diag_pre_multiply(sigma_G, L_rho) * z[j])[1];
    bMY_G[j] = bMY + (diag_pre_multiply(sigma_G, L_rho) * z[j])[2];
    
    } 
}
model{
    vector[360] mu;
    L_rho ~ lkj_corr_cholesky( 2 );
    sigma ~ exponential( 1 );
    sigma_G ~ exponential( 1 );
    bMY ~ normal( -1 , 0.5 );
    a ~ normal( 5 , 2 );
    
    for (i in 1:2) {
    z[i] ~ normal(0, 1);  // This applies the normal distribution to each 2D vector
    }

    for ( i in 1:360 ) {
        mu[i] = a_G[G[i]] + bMY_G[G[i]] * MY_std[i];
    }
    GS_std ~ normal( mu , sigma );
}

generated quantities{
    vector[360] log_lik;
    vector[360] mu;
    matrix[2,2] Rho;
    Rho = multiply_lower_tri_self_transpose(L_rho);
    for ( i in 1:360 ) {
        mu[i] = a_G[G[i]] + bMY_G[G[i]] * MY_std[i];
    }
    for ( i in 1:360 ) log_lik[i] = normal_lpdf( GS_std[i] | mu[i] , sigma );
}

