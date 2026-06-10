data {
  int<lower=1> N;                       
  int<lower=1> S;   
  array[N] int subject; 
  array[N] int trial; // 1 if first trial of block
  array[N] real b1;                       
  array[N] real b2;                      
  array[N] int<lower=0,upper=1> choice;   
}
 
parameters {
  // Independent parameters for each subject on unconstrained space
  vector[S] logit_alpha_pos;
  vector[S] logit_alpha_neg;
  vector[S] log_beta;
}

transformed parameters {
  // Transform to natural scale
  vector<lower=0,upper=1>[S] alpha_pos = inv_logit(logit_alpha_pos);
  vector<lower=0,upper=1>[S] alpha_neg = inv_logit(logit_alpha_neg);
  vector<lower=0>[S] beta              = exp(log_beta);
}

model {
  // Weakly informative priors for each individual subject
  // These replace the hierarchical group distributions
  logit_alpha_pos ~ normal(0, 1.5); 
  logit_alpha_neg ~ normal(0, 1.5);
  log_beta        ~ normal(0, 2);
  
  // initialize Q-values
  matrix[2, N+1] Q = rep_matrix(0, 2, N+1);
  vector[N] logits;

  for (i in 1:N){
    int s = subject[i]; // current subject index

    if (trial[i] == 1){ 
      Q[,i] = rep_vector(0.5, 2); // reset first trial of block
    }

    // choice probability logic
    logits[i] = beta[s] * (Q[1,i] - Q[2,i]);

    // prediction errors
    real d1 = b1[i] - Q[1,i];
    real d2 = b2[i] - Q[2,i];

    // learning rate selection
    real lr1 = d1 > 0 ? alpha_pos[s] : alpha_neg[s];
    real lr2 = d2 > 0 ? alpha_pos[s] : alpha_neg[s];

    // update Q-values
    Q[1,i+1] = Q[1,i] + lr1 * d1;
    Q[2,i+1] = Q[2,i] + lr2 * d2;
  }

  // likelihood
  choice ~ bernoulli_logit(logits);
} 

generated quantities {
  vector[S] alpha_base = 0.5 * (alpha_pos + alpha_neg);
  vector[S] bias       = alpha_pos ./ (alpha_pos + alpha_neg);
}
 