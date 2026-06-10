rm(list = ls())
gc()

{
  source("Model Fitting/setup.R")
  data_dir = file.path("Data Simulation", "Simulated Data")
  fit_dir = file.path("Model Fitting", "Fit Models")
  save_dir =  file.path("Simulations and Analysis", "RDS files")
  if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
} # setup

{
  simulate_const = function(fit, stan_data, b1_better){
    
    d = stan_data
    
    {
      N = d$N
      S = d$S
      subject = d$subject
      b1 = d$b1
      b2 = d$b2
      choice = d$choice
      trial = d$trial
      
    } # unpack data
    {
      draws <- as_draws_df(fit)
      D <- nrow(draws)           # number of posterior draws
      
      beta  <- as.matrix(fit, pars = "beta")
      alpha <- as.matrix(fit, pars = "alpha") 
      
    } # unpack draws
    {
      Q1 <- matrix(0.5, D, N+1)
      Q2 <- matrix(0.5, D, N+1)
      
      acc_accumulator <- matrix(0, D, max(trial))
      stay_r_sum  <- rep(0, D)
      stay_l_sum  <- rep(0, D)
      r_count     <- rep(0, D)
      l_count     <- rep(0, D)
      
      pb <- progress_bar$new(
        total = N,
        format = "  Constant [1/3] [:bar] :percent eta: :eta",
        stream = stderr() 
      )
    } # initialize
    for(i in 1:N){
      pb$tick()
      
      
      if(trial[i] == 1){
        Q1[,i] <- 0.5
        Q2[,i] <- 0.5
        
      }
      
      s <- subject[i]
      
      logit <- beta[,s] * (Q1[,i] - Q2[,i])
      
      p <- 1/(1 + exp(-logit))
      choice_sim <- rbinom(D, 1, p)
      
      is_correct <- (choice_sim == b1_better[i])
      acc_accumulator[, trial[i]] <- acc_accumulator[, trial[i]] + is_correct
      
      # stay / switch calculation
      {if(trial[i] > 1){
        stay <- as.numeric(choice_sim == prev_choice)
        # reward
        idx_r <- which(prev_r == 1)
        stay_r_sum[idx_r] <- stay_r_sum[idx_r] + stay[idx_r]
        r_count[idx_r] <- r_count[idx_r] + 1
        
        # loss 
        idx_l <- which(prev_r == 0)
        stay_l_sum[idx_l] <- stay_l_sum[idx_l] + stay[idx_l]
        l_count[idx_l] <- l_count[idx_l] + 1
      }
        prev_choice = choice_sim
        prev_r = ifelse(choice_sim==1, b1[i], b2[i])
      }
      
      d1 <- ifelse(choice[i] == 1,  b1[i] - Q1[,i] , 0)
      d2 <- ifelse(choice[i] == 0,  b2[i] - Q2[,i] , 0) 
      
      Q1[,i+1] <- Q1[,i] + alpha[,s] * d1
      Q2[,i+1] <- Q2[,i] + alpha[,s] * d2
      
      
    }  # learning loop
    
    acc_per_trial_draw <- acc_accumulator / (N/max(trial))
    
    stay_p = data.frame(after_reward = stay_r_sum / r_count,
                        after_loss   = stay_l_sum / l_count,
                        general      = (stay_r_sum + stay_l_sum) / (r_count + l_count))
    
    
    return(list(acc_per_trial_draw,stay_p))
  }
  simulate_asymm = function(fit, stan_data, b1_better){
    
    d = stan_data
    
    {
      N = d$N
      S = d$S
      subject = d$subject
      b1 = d$b1
      b2 = d$b2
      choice = d$choice
      trial = d$trial
    } # unpack data
    {
      draws <- as_draws_df(fit)
      D <- nrow(draws)           # number of posterior draws
      
      beta  <- as.matrix(fit, pars = "beta")
      alpha_pos <- as.matrix(fit, pars = "alpha_pos")
      alpha_neg <- as.matrix(fit, pars = "alpha_neg")
      
    } # unpack draws
    {
      Q1 <- matrix(0.5, D, N+1)
      Q2 <- matrix(0.5, D, N+1)
      
      acc_accumulator <- matrix(0, D, max(trial))
      
      stay_r_sum  <- rep(0, D)
      stay_l_sum  <- rep(0, D)
      r_count     <- rep(0, D)
      l_count     <- rep(0, D)
      
      pb <- progress_bar$new(
        total = N,
        format = "  Asymmetric [2/3] [:bar] :percent eta: :eta"
      )
    } # initialize
    for(i in 1:N){
      pb$tick()
      
      
      if(trial[i] == 1){
        Q1[,i] <- 0.5
        Q2[,i] <- 0.5
        
      }
      
      s <- subject[i]
      
      logit <- beta[,s] * (Q1[,i] - Q2[,i])
      
      p <- 1/(1 + exp(-logit))
      choice_sim <- rbinom(D, 1, p)
      
      is_correct <- (choice_sim == b1_better[i])
      acc_accumulator[, trial[i]] <- acc_accumulator[, trial[i]] + is_correct
      
      # stay / switch calculation
      {if(trial[i] > 1){
        stay <- as.numeric(choice_sim == prev_choice)
        # reward
        idx_r <- which(prev_r == 1)
        stay_r_sum[idx_r] <- stay_r_sum[idx_r] + stay[idx_r]
        r_count[idx_r] <- r_count[idx_r] + 1
        
        # loss 
        idx_l <- which(prev_r == 0)
        stay_l_sum[idx_l] <- stay_l_sum[idx_l] + stay[idx_l]
        l_count[idx_l] <- l_count[idx_l] + 1
      }
        prev_choice = choice_sim
        prev_r = ifelse(choice_sim==1, b1[i], b2[i])
      }
      
      d1 <- ifelse(choice[i] == 1,  b1[i] - Q1[,i] , 0)
      d2 <- ifelse(choice[i] == 0,  b2[i] - Q2[,i] , 0) 
      
      lr1 <- ifelse(d1 > 0, alpha_pos[,s], alpha_neg[,s])
      lr2 <- ifelse(d2 > 0, alpha_pos[,s], alpha_neg[,s])
      
      Q1[,i+1] <- Q1[,i] + lr1 * d1
      Q2[,i+1] <- Q2[,i] + lr2 * d2
      
      
    }  # learning loop
    
    acc_per_trial_draw <- acc_accumulator / (N/max(trial))
    
    stay_p = data.frame(after_reward = stay_r_sum / r_count,
                        after_loss   = stay_l_sum / l_count,
                        general      = (stay_r_sum + stay_l_sum) / (r_count + l_count))
    
    
    return(list(acc_per_trial_draw,stay_p))
  }
  simulate_true = function(fit, stan_data, b1_better){
    
    d = stan_data
    
    {
      N = d$N
      S = d$S
      subject = d$subject
      b1 = d$b1
      b2 = d$b2
      choice = d$choice
      trial = d$trial
    } # unpack data
    {
      draws <- as_draws_df(fit)
      D <- nrow(draws)           # number of posterior draws
      
      beta  <- as.matrix(fit, pars = "beta")
    } # unpack draws
    {
      Q1 <- matrix(0.5, D, N+1)
      Q2 <- matrix(0.5, D, N+1)
      
      acc_accumulator <- matrix(0, D, max(trial))
      stay_r_sum  <- rep(0, D)
      stay_l_sum  <- rep(0, D)
      r_count     <- rep(0, D)
      l_count     <- rep(0, D)
      
      pb <- progress_bar$new(
        total = N,
        format = "  Diminishing [3/3] [:bar] :percent eta: :eta"
      )
    } # initialize
    for(i in 1:N){
      pb$tick()
      
      
      if(trial[i] == 1){
        Q1[,i] <- 0.5
        Q2[,i] <- 0.5
        
        counts1 <- 1
        counts2 <- 1
      }
      
      s <- subject[i]
      
      logit <- beta[,s] * (Q1[,i] - Q2[,i])
      
      p <- 1/(1 + exp(-logit))
      choice_sim <- rbinom(D, 1, p) # 0: b2, 1: b1
      
      is_correct <- (choice_sim == b1_better[i])
      acc_accumulator[, trial[i]] <- acc_accumulator[, trial[i]] + is_correct
      
      # stay / switch calculation
      {if(trial[i] > 1){
        stay <- as.numeric(choice_sim == prev_choice)
        # reward
        idx_r <- which(prev_r == 1)
        stay_r_sum[idx_r] <- stay_r_sum[idx_r] + stay[idx_r]
        r_count[idx_r] <- r_count[idx_r] + 1
        
        # loss 
        idx_l <- which(prev_r == 0)
        stay_l_sum[idx_l] <- stay_l_sum[idx_l] + stay[idx_l]
        l_count[idx_l] <- l_count[idx_l] + 1
        }
      prev_choice = choice_sim
      prev_r = ifelse(choice_sim==1, b1[i], b2[i])
      }
      
      d1 <- ifelse(choice[i] == 1,  b1[i] - Q1[,i] , 0)
      d2 <- ifelse(choice[i] == 0,  b2[i] - Q2[,i] , 0) 
      
      Q1[,i+1] <- Q1[,i] + 1.0/counts1 * d1
      Q2[,i+1] <- Q2[,i] + 1.0/counts2 * d2
      
      
      counts1 = counts1 + ifelse(choice[i] == 1, 1 , 0)
      counts2 = counts2 + ifelse(choice[i] == 0, 1 , 0) 
      
      
    }  # learning loop
    
    acc_per_trial_draw <- acc_accumulator / (N/max(trial))
    
    stay_p = data.frame(after_reward = stay_r_sum / r_count,
                        after_loss   = stay_l_sum / l_count,
                        general      = (stay_r_sum + stay_l_sum) / (r_count + l_count))
    
    
    return(list(acc_per_trial_draw,stay_p))
  }
  
  summarize_acc <- function(acc_matrix, label) {
    data.frame(
      trial = 1:ncol(acc_matrix),
      mean  = colMeans(acc_matrix),
      lo    = apply(acc_matrix, 2, quantile, 0.05),
      hi    = apply(acc_matrix, 2, quantile, 0.95),
      model = label
    )
  }
} # functions

suppressMessages({
  data =  read_csv(file.path(data_dir,"Case_2.csv"))
  
  stan_data = fetch_stan_data("Case_2")
  
  block_summ = data %>%
    group_by(subject, block) %>%
    summarise(b1 = mean(b1),
              b2 = mean(b2)) %>%
    mutate(b1_better = b1 > b2)
  
  data = data %>%
    left_join(block_summ %>% select(subject,block, b1_better)) %>%
    mutate(correct = as.integer(b1_better) == choice ) 
  
  trial_summ = data %>%
    group_by(trial) %>%
    summarise(accuracy = mean(correct))
  
  b1_better = data$b1_better
}) # import and summarize data 

{
  fit_const = readRDS(file.path(fit_dir, "Case_2_const.rds"))
  sim_const = simulate_const(fit_const, stan_data, b1_better) # [[1]] accuracy df  [[2]] stay_p df
  
  fit_asymm = readRDS(file.path(fit_dir, "Case_2_asymm.rds"))
  sim_asymm = simulate_asymm(fit_asymm, stan_data, b1_better)
  
  fit_true = readRDS(file.path(fit_dir, "Case_2_true.rds"))
  sim_true = simulate_true(fit_true, stan_data, b1_better)
} # Simulate

{
  ppc_const <- summarize_acc(sim_const[[1]], "Constant LR")
  ppc_asymm   <- summarize_acc(sim_asymm[[1]],  "Asymmetric LR")
  ppc_true  <- summarize_acc(sim_true[[1]], "True")
  
  acc_ppc_case_2 <- bind_rows(ppc_const, ppc_asymm, ppc_true)
} # Calculate accuracy PPC

{
stay_ppc_case_2 =  bind_rows(sim_const[[2]] %>% as.data.frame() %>% mutate(model = "Constant LR"),
                             sim_asymm[[2]] %>% as.data.frame() %>% mutate(model = "Asymmetric LR"),
                             sim_true[[2]]  %>% as.data.frame() %>% mutate(model = "True"))

} # calculate stay PPC

# Save PPC results
saveRDS(acc_ppc_case_2, file =  file.path(save_dir,"acc_ppc_case_2.rds") )
saveRDS(stay_ppc_case_2, file =  file.path(save_dir,"stay_ppc_case_2.rds") )


if (F) {
  
{
  ggplot(acc_ppc_case_2, aes(x = trial, y = mean, color = model, fill = model)) +
    # Credible interval ribbons
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
    # Mean lines
    geom_line(linewidth = 1) +
    # Real data (red points)
    geom_point(data = trial_summ, aes(x = trial, y = accuracy), 
               inherit.aes = FALSE, color = "black", size = 2) +
    theme_bw() +
    labs(title = "Model Comparison: Accuracy per trial",
         subtitle = "Blue/Green = Models, Red = Observed Data",
         y = "Mean Accuracy",
         x = "trial")
} # plot summary acc

{
  stay_ppc = data %>%
    group_by(block) %>%
    mutate(stay = choice == lag(choice)) %>%
    ungroup() %>%
    summarise(stay_p = mean(stay,na.rm = T))

  
  ggplot(stay_ppc_case_2)+ 
    geom_vline(xintercept = stay_ppc$stay_p, linetype = 2)+
    geom_histogram(aes(x = general, fill = model),alpha = 0.5, position = "identity",)+
    theme_bw()

  } # plot stay - switch probabilities

  } # Visualize PPC results