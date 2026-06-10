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
  // group-level location and scale on unconstrained space
  vector[3] mu;                           // means for [logit(alpha_pos), logit(alpha_neg), log(beta)]
  vector<lower=1e-6>[3] sigma;            // scales
  cholesky_factor_corr[3] L_Omega;        // Cholesky factor of correlation

  // subject effects (standard normal)
  matrix[3, S] z;
}

transformed parameters {
  // subject parameters on unconstrained space
  matrix[3, S] theta;
  theta = diag_pre_multiply(sigma, L_Omega) * z + rep_matrix(mu, S);
  
  // subject parameters on natural scale
  vector[S] alpha_pos = inv_logit(theta[1, ]');  
  vector[S] alpha_neg = inv_logit(theta[2, ]');  
  vector[S] beta      = exp(theta[3, ]');  
}


model {
  // priors
  mu ~ normal(0, 1);          
  sigma ~ lognormal(0, 1);    
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(z) ~ normal(0, 1);
  
  // initialize Q-values
  matrix[2, N+1] Q = rep_matrix(0, 2, N+1);
  vector[N] logits;

  for (i in 1:N){
    if (trial[i] == 1){ 
      Q[,i] = rep_vector(0.5, 2); // first trial of block
    }

    // choice probability
    logits[i] = beta[subject[i]] * (Q[1,i] - Q[2,i]);

    // prediction errors
    real d1 = b1[i] - Q[1,i];
    real d2 = b2[i] - Q[2,i];

    // select learning rate based on sign of PE
    real lr1 = d1 > 0 ? alpha_pos[subject[i]] : alpha_neg[subject[i]];
    real lr2 = d2 > 0 ? alpha_pos[subject[i]] : alpha_neg[subject[i]];

    // update Q-values
    Q[1,i+1] = Q[1,i] + lr1 * d1;
    Q[2,i+1] = Q[2,i] + lr2 * d2;
  }

  // likelihood
  choice ~ bernoulli_logit(logits);
} 

generated quantities {
  // recover correlation matrix
  corr_matrix[3] Omega = multiply_lower_tri_self_transpose(L_Omega);

  vector[S] alpha_base = 0.5 * (alpha_pos + alpha_neg);
  vector[S] bias       = alpha_pos ./ (alpha_pos + alpha_neg);
}
