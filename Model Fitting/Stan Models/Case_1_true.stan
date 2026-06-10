data {
  int<lower=1> N;                       
  int<lower=1> S;   
  array[N] int subject; 
  array[N] int trial; 
  array[N] real b1;                       
  array[N] real b2;                      
  array[N] int<lower=0,upper=1> choice;   
}
 
parameters {
  // group level location and scale on unconstrained space
  vector[1] mu;                           // means for [logit(alpha), log(beta)]
  vector<lower=1e-6>[1] sigma;            // scales
  cholesky_factor_corr[1] L_Omega;        // Cholesky factor of correlation

  // subject effects (standard normal)
  matrix[1, S] z;
}

transformed parameters {
  // subject parameters on unconstrained space
  matrix[1, S] theta;
  theta = diag_pre_multiply(sigma, L_Omega) * z + rep_matrix(mu, S);
  
  vector[S] beta  = exp(theta[1, ]');        // theta[2, ] is row for beta
}

model {
  // priors
  mu ~ normal(0, 1);          // weakly informative
  sigma ~ lognormal(0, 1);    
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(z) ~ normal(0, 1);
  
  matrix[2, N+1] Q;
  vector[N] logits;
  real lr;

  
  Q = rep_matrix(0, 2, N+1);  // 2 rows, N columns of zeros

  for (i in 1:N){
    
    if (trial[i]==1){ 
        Q[,i] = rep_vector(0.5, 2); // set Q values to 0.5 for first trials
    }
    
    logits[i] = beta[subject[i]] * trial[i] * (Q[1,i] - Q[2,i]);
    
    lr = 1.0 / (trial[i]);
    
    Q[1,i+1] = Q[1,i] + lr  * (b1[i] - Q[1,i]);
    Q[2,i+1] = Q[2,i] + lr  * (b2[i] - Q[2,i]);
  }
  
  choice ~ bernoulli_logit(logits);
} 

generated quantities {
  // recover correlation and subject-level parameters on natural scales
  corr_matrix[1] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
