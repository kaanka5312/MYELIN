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
    vector[2] z[20];
}

//transformed parameters{
//    vector[20] b_cafe;
//    vector[20] a_cafe;
//    a_cafe = (diag_pre_multiply(sigma_cafe, L_rho) * z_alfa);
//    b_cafe = (diag_pre_multiply(sigma_cafe, L_rho) * z_beta);
    //}

transformed parameters{
    vector[20] b_cafe;
    vector[20] a_cafe;
    for ( j in 1:20 ){
    a_cafe[j] = a + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[1];
    b_cafe[j] = b + (diag_pre_multiply(sigma_cafe, L_rho) * z[j])[2];
    }
}

model{
    vector[200] mu;
    L_rho ~ lkj_corr_cholesky( 2 );
    sigma ~ exponential( 1 );
    sigma_cafe ~ exponential( 1 );
    b ~ normal( -1 , 0.5 );
    a ~ normal( 5 , 2 );
    // Properly specify prior for each element of each vector in the z array
    for (i in 1:20) {
        z[i] ~ normal(0, 1);  // This applies the normal distribution to each 2D vector
    }
//    to_vector( z ) ~ normal( 0 , 1 );
//  {
//    vector[2] YY[20];
//    vector[2] MU;
//    MU = [ a , b ]';
//    for ( j in 1:20 ) YY[j] = [ a_cafe[j] , b_cafe[j] ]';
//    YY ~ multi_normal( MU , quad_form_diag(Rho , sigma_cafe) );
//    }
    
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

