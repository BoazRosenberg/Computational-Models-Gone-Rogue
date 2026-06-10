rm(list = ls())
gc()

{
  source("Model Fitting/setup.R")
  data_dir = file.path("Data Simulation", "Simulated Data")
  fit_dir = file.path("Model Fitting", "Fit Models")
  save_dir =  file.path("Simulations and Analysis", "RDS files")
  if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
} # setup

##### Simulate #####

{
  log_mean_exp <- function(x) {
    max_x <- max(x)
    return(max_x + log(mean(exp(x - max_x))))
  }
  
  simulate_const = function(fit, stan_data){
    
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
      
      alpha <- as.matrix(fit, pars = "alpha")
      beta  <- as.matrix(fit, pars = "beta")
      
    } # unpack draws
    {
      Q1 <- matrix(0.5, D, N+1)
      Q2 <- matrix(0.5, D, N+1)
      
      loglik <- matrix(0, D, N)
      
      PE     <-  matrix(NA, D, N) 
      update <-  matrix(NA, D, N)
      
      pb <- progress_bar$new(
        total = N,
        format = "  Constant [1/3] [:bar] :percent eta: :eta",
        stream = stderr(),
        force = TRUE
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
      
      loglik[,i] <- dbinom(choice[i], 1, p, log = TRUE)
      
      d1 <- ifelse(choice[i] == 1,  b1[i] - Q1[,i] , 0)
      d2 <- ifelse(choice[i] == 0,  b2[i] - Q2[,i] , 0) 
      
      Q1[,i+1] <- Q1[,i] + alpha[,s] * d1
      Q2[,i+1] <- Q2[,i] + alpha[,s] * d2
      
      # record PEs and updates
      PE[,i]     <-  ifelse(choice[i]==1, d1, d2)
      update[,i] <-  ifelse(choice[i]==1, alpha[,s] * d1, alpha[,s] * d2)
      
    }  # learning loop
    
    # Per round means
    {
      # Calculate the Log Predictive Density (LPD) per trial
      lpd_trial <- apply(loglik, 2, log_mean_exp)   
      # Average the trial-wise LPDs by round
      mean_lpd_round <- tapply(lpd_trial, trial, mean)
      
      #mean_loglik_trial <- colMeans(loglik)
      #mean_loglik_round <- tapply(mean_loglik_trial, trial, mean)
    } # mean likelihood
    
    {
      mean_PE_trial <- colMeans(PE, na.rm = TRUE) # over draws
      mean_PE_round <- tapply(mean_PE_trial, trial, mean) # over subjects and blocks
      
      # Separate by positive and negative
      PE_pos <- ifelse(PE > 0, PE, NA)
      PE_neg <- ifelse(PE < 0, PE, NA)
      
      mean_PE_pos_trial <- colMeans(PE_pos, na.rm = TRUE)
      mean_PE_neg_trial <- colMeans(PE_neg, na.rm = TRUE)
      
      mean_PE_pos_round <- tapply(mean_PE_pos_trial, trial, mean, na.rm = TRUE)
      mean_PE_neg_round <- tapply(mean_PE_neg_trial, trial, mean, na.rm = TRUE)
      
    } # mean prediction error
    {
      mean_update_trial <- colMeans(update, na.rm = TRUE) # over draws
      mean_update_round <- tapply(mean_update_trial, trial, mean) # over subjects and blocks
      
      # Separate by positive and negative
      update_pos <- ifelse(PE > 0, update, NA)
      update_neg <- ifelse(PE < 0, update, NA)
      
      mean_update_pos_trial <- colMeans(update_pos, na.rm = TRUE)
      mean_update_neg_trial <- colMeans(update_neg, na.rm = TRUE)
      
      mean_update_pos_round <- tapply(mean_update_pos_trial, trial, mean, na.rm = TRUE)
      mean_update_neg_round <- tapply(mean_update_neg_trial, trial, mean, na.rm = TRUE)
    } # mean update
    
    
    df = data.frame(trial = 1:30,
                    ll = mean_lpd_round,
                    PE = mean_PE_round,
                    update = mean_update_round,
                    
                    PE_pos = mean_PE_pos_round,
                    PE_neg = mean_PE_neg_round,
                    
                    update_pos = mean_update_pos_round,
                    update_neg = mean_update_neg_round)
    
    return(df)
  }
  
  simulate_asymm = function(fit, stan_data){
    
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
      
      alpha_pos <- as.matrix(fit, pars = "alpha_pos")
      alpha_neg <- as.matrix(fit, pars = "alpha_neg")
      beta      <- as.matrix(fit, pars = "beta")
    } # unpack draws
    {
      Q1 <- matrix(0.5, D, N+1)
      Q2 <- matrix(0.5, D, N+1)
      
      loglik <- matrix(0, D, N)
      
      PE     <-  matrix(NA, D, N) 
      update <-  matrix(NA, D, N)
      
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
      
      logit <-  beta[,s] * (Q1[,i] - Q2[,i])
      
      p <- 1/(1 + exp(-logit))
      
      loglik[,i] <- dbinom(choice[i], 1, p, log = TRUE)
      
      d1 <- ifelse(choice[i] == 1,  b1[i] - Q1[,i] , 0)
      d2 <- ifelse(choice[i] == 0,  b2[i] - Q2[,i] , 0) 
      
      lr1 <- ifelse(d1 > 0, alpha_pos[,s], alpha_neg[,s])
      lr2 <- ifelse(d2 > 0, alpha_pos[,s], alpha_neg[,s])
      
      Q1[,i+1] <- Q1[,i] + lr1 * d1
      Q2[,i+1] <- Q2[,i] + lr2 * d2
      
      # record PEs and updates
      PE[,i]     <-  ifelse(choice[i]==1,
                            d1,
                            d2)
      update[,i] <-  ifelse(choice[i]==1,
                            lr1 * d1,
                            lr2 * d2)
      
      
    }  # learning loop
    
    {
      # Calculate the Log Predictive Density (LPD) per trial
      lpd_trial <- apply(loglik, 2, log_mean_exp)   
      # Average the trial-wise LPDs by round
      mean_lpd_round <- tapply(lpd_trial, trial, mean)
      
      #mean_loglik_trial <- colMeans(loglik)
      #mean_loglik_round <- tapply(mean_loglik_trial, trial, mean)
    } # mean likelihood
    {
      mean_PE_trial <- colMeans(PE, na.rm = TRUE) # over draws
      mean_PE_round <- tapply(mean_PE_trial, trial, mean) # over subjects and blocks
      
      # Separate by positive and negative
      PE_pos <- ifelse(PE > 0, PE, NA)
      PE_neg <- ifelse(PE < 0, PE, NA)
      
      mean_PE_pos_trial <- colMeans(PE_pos, na.rm = TRUE)
      mean_PE_neg_trial <- colMeans(PE_neg, na.rm = TRUE)
      
      mean_PE_pos_round <- tapply(mean_PE_pos_trial, trial, mean, na.rm = TRUE)
      mean_PE_neg_round <- tapply(mean_PE_neg_trial, trial, mean, na.rm = TRUE)
      
    } # mean prediction error
    {
      mean_update_trial <- colMeans(update, na.rm = TRUE) # over draws
      mean_update_round <- tapply(mean_update_trial, trial, mean) # over subjects and blocks
      
      # Separate by positive and negative
      update_pos <- ifelse(PE > 0, update, NA)
      update_neg <- ifelse(PE < 0, update, NA)
      
      mean_update_pos_trial <- colMeans(update_pos, na.rm = TRUE)
      mean_update_neg_trial <- colMeans(update_neg, na.rm = TRUE)
      
      mean_update_pos_round <- tapply(mean_update_pos_trial, trial, mean, na.rm = TRUE)
      mean_update_neg_round <- tapply(mean_update_neg_trial, trial, mean, na.rm = TRUE)
    } # mean update
    
    df = data.frame(trial = 1:30,
                    ll = mean_lpd_round,
                    PE = mean_PE_round,
                    update = mean_update_round,
                    
                    PE_pos = mean_PE_pos_round,
                    PE_neg = mean_PE_neg_round,
                    
                    update_pos = mean_update_pos_round,
                    update_neg = mean_update_neg_round)
    
    return(df)
  }
  
  simulate_dim = function(fit, stan_data){
    
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
      
      counts1 <- 1
      counts2 <- 1
      
      loglik <- matrix(0, D, N)
      
      PE     <-  matrix(NA, D, N) 
      update <-  matrix(NA, D, N)
      
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
      
      loglik[,i] <- dbinom(choice[i], 1, p, log = TRUE)
      
      d1 <- ifelse(choice[i] == 1,  b1[i] - Q1[,i] , 0)
      d2 <- ifelse(choice[i] == 0,  b2[i] - Q2[,i] , 0) 
      
      Q1[,i+1] <- Q1[,i] + 1/counts1 * d1
      Q2[,i+1] <- Q2[,i] + 1/counts2* d2
      
      
      # record PEs and updates
      PE[,i]     <-  ifelse(choice[i]==1, d1, d2)
      update[,i] <-  ifelse(choice[i]==1, 
                            1/counts1 * d1,
                            1/counts2 * d2)
      
      
      counts1 = counts1 + ifelse(choice[i] == 1 ,1, 0)
      counts2 = counts2 + ifelse(choice[i] == 0 ,1, 0)
      
    }  # learning loop
    
    # Per round means
    {
      # Calculate the Log Predictive Density (LPD) per trial
      lpd_trial <- apply(loglik, 2, log_mean_exp)   
      # Average the trial-wise LPDs by round
      mean_lpd_round <- tapply(lpd_trial, trial, mean)
      
      #mean_loglik_trial <- colMeans(loglik)
      #mean_loglik_round <- tapply(mean_loglik_trial, trial, mean)
    } # mean likelihood
    
    {
      mean_PE_trial <- colMeans(PE, na.rm = TRUE) # over draws
      mean_PE_round <- tapply(mean_PE_trial, trial, mean) # over subjects and blocks
      
      # Separate by positive and negative
      PE_pos <- ifelse(PE > 0, PE, NA)
      PE_neg <- ifelse(PE < 0, PE, NA)
      
      mean_PE_pos_trial <- colMeans(PE_pos, na.rm = TRUE)
      mean_PE_neg_trial <- colMeans(PE_neg, na.rm = TRUE)
      
      mean_PE_pos_round <- tapply(mean_PE_pos_trial, trial, mean, na.rm = TRUE)
      mean_PE_neg_round <- tapply(mean_PE_neg_trial, trial, mean, na.rm = TRUE)
      
    } # mean prediction error
    {
      mean_update_trial <- colMeans(update, na.rm = TRUE) # over draws
      mean_update_round <- tapply(mean_update_trial, trial, mean) # over subjects and blocks
      
      # Separate by positive and negative
      update_pos <- ifelse(PE > 0, update, NA)
      update_neg <- ifelse(PE < 0, update, NA)
      
      mean_update_pos_trial <- colMeans(update_pos, na.rm = TRUE)
      mean_update_neg_trial <- colMeans(update_neg, na.rm = TRUE)
      
      mean_update_pos_round <- tapply(mean_update_pos_trial, trial, mean, na.rm = TRUE)
      mean_update_neg_round <- tapply(mean_update_neg_trial, trial, mean, na.rm = TRUE)
    } # mean update
    
    
    df = data.frame(trial = 1:30,
                    ll = mean_lpd_round,
                    PE = mean_PE_round,
                    update = mean_update_round,
                    
                    PE_pos = mean_PE_pos_round,
                    PE_neg = mean_PE_neg_round,
                    
                    update_pos = mean_update_pos_round,
                    update_neg = mean_update_neg_round)
    
    return(df)
  }
  
} # functions

stan_data_partial = fetch_stan_data("Case_2")

{
  df_asymm = simulate_asymm(fit = readRDS(file.path(fit_dir, "Case_2_asymm.rds")),
                            stan_data = stan_data_partial)
  
  df_const = simulate_const(fit = readRDS(file.path(fit_dir, "Case_2_const.rds")),
                            stan_data = stan_data_partial)
  
  df_true = simulate_dim(fit = readRDS(file.path(fit_dir, "Case_2_true.rds")),
                        stan_data = stan_data_partial)
  
} # simulate


update_simulation_case_2 = bind_rows(df_asymm %>% mutate(model = "asymm"),
                                     df_const%>% mutate(model = "const"),
                                     df_true %>% mutate(model = "true"))


saveRDS(update_simulation_case_2, file = file.path(save_dir,"sim_case_2.rds"))


##### Visualize ####

if (F){
means = update_simulation_case_2 %>% group_by(model) %>% summarise(mean_ll = mean(ll))

{ggplot(update_simulation_case_2, aes(x = trial, col = model))+
  geom_point(aes(y = exp(ll)))+
  geom_smooth(aes(y = exp(ll)), se = F)+
  geom_hline(data = means, aes(yintercept = exp(mean_ll), col = model))+
  theme_bw()+
  labs(x = "Trial", y = "Average Likelihood")} # Likelihood

{ggplot(update_simulation_case_2 %>% 
         select(model, trial, PE_pos, PE_neg) %>% 
         gather(key = "direction", value = "PE", PE_pos, PE_neg)%>%
          mutate(direction = factor(direction, levels = c("PE_pos", "PE_neg"))),
       aes(x = trial, col = model))+
  geom_line(aes(y = abs(PE))) +
  facet_grid(direction~.)+
  theme_bw()+
  labs(x = "Trial")} # PE
  
{ggplot(update_simulation_case_2 %>% 
          select(model, trial, update_pos, update_neg) %>% 
          gather(key = "direction", value = "update", update_pos, update_neg) %>%
          mutate(direction = factor(direction, levels = c("update_pos", "update_neg"))),
        aes(x = trial, col = model))+
    geom_smooth(aes(y = abs(update)), se = F) +
    facet_grid(direction~.)+
    theme_bw()+
    labs(x = "Trial")} # update
}
