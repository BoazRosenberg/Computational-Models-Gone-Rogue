suppressMessages({
library(readr)
library(posterior)
library(dplyr)
library(tidyr)
library(shinystan)
library(rstan)
library(ggplot2)
library(progress)
library(parallel)
})

options(dplyr.summarise.inform = FALSE)

fetch_stan_data = function(model_name){
  
  dir = file.path("Data Simulation",
                  "Simulated Data")
  
  df = suppressMessages({read_csv(paste0(dir,"/",model_name,".csv"))})
  
  
  stan_data <- list(
    N = nrow(df),
    S = length(unique(df$subject)),
    subject = as.integer(factor(df$subject)),
    trial       = as.integer(factor(df$trial)),
    choice  = df$choice,
    b1 = df$b1,
    b2 = df$b2
  )
  
  return(stan_data)
  
}

fit_and_save <- function(name, stan_data, iterations = 2500, warmup = 1000, suffix = "") {
  
  # Compile the model
  mod <- stan_model(file = file.path("Model Fitting",
                                     "Stan Models",
                                     paste0(name, ".stan")),
                    model_name = name)
  
  # Sample from the posterior
  fit <- sampling(
    object = mod,
    data = stan_data,
    chains = 4,
    cores = detectCores(),
    iter = iterations,             
    warmup = warmup,
    show_messages = TRUE,
    
  )
  
  # Save the fit and the compiled model
  
  dir = file.path("Model Fitting",
                  "Fit Models")
  
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  file_name = file.path(dir, paste0(name, suffix, ".rds"))
  
  saveRDS(fit, file = file_name)
  
  return(fit)
}

