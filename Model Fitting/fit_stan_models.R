
source("Model Fitting/setup.R")

{
  Case_1_stan_data = fetch_stan_data("Case_1")
  Case_2_stan_data = fetch_stan_data("Case_2")
  
  Case_2_extra_opposing_data = fetch_stan_data("Case_2_extra_opposing")
  Case_2_extra_random_data = fetch_stan_data("Case_2_extra_random")

  
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)  # caches compiled models
  
  {
    log_progress <- function(msg) {
      now <- Sys.time()
      elapsed <- as.numeric(difftime(now, start_time, units = "mins"))
      cat(sprintf("[%s | +%.1f min] %s\n", format(now, "%H:%M"), elapsed, msg))
    }
  } # log progress function

} # setup

{
  log_progress("[1/9] Fitting Case_1_const")
  Case_1_const = fit_and_save("Case_1_const", Case_1_stan_data)
  Case_1_const_bridge = bridgesampling::bridge_sampler(Case_1_const, silent = TRUE)
  saveRDS(Case_1_const_bridge,
          file = file.path("Model Fitting", "Fit Models","Case_1_const_bridge.rds"))
  
  log_progress("[2/9] Fitting Case_1_dim")
  Case_1_dim = fit_and_save("Case_1_dim", Case_1_stan_data)
  Case_1_dim_bridge = bridgesampling::bridge_sampler(Case_1_dim, silent = TRUE)
  saveRDS(Case_1_dim_bridge,
          file = file.path("Model Fitting", "Fit Models","Case_1_dim_bridge.rds"))
  
  log_progress("[3/9] Fitting Case_1_true")
  Case_1_true = fit_and_save("Case_1_true", Case_1_stan_data)
  Case_1_true_bridge = bridgesampling::bridge_sampler(Case_1_true, silent = TRUE)
  saveRDS(Case_1_true_bridge,
          file = file.path("Model Fitting", "Fit Models","Case_1_true_bridge.rds"))
  
} # case 1

{
  log_progress("[4/9] Fitting Case_2_asymm")
  Case_2_asymm = fit_and_save("Case_2_asymm", Case_2_stan_data)
  Case_2_asymm_bridge = bridgesampling::bridge_sampler(Case_2_asymm, silent = TRUE)
  saveRDS(Case_2_asymm_bridge,
          file = file.path("Model Fitting", "Fit Models","Case_2_asymm_bridge.rds"))
  
  log_progress("[5/9] Fitting Case_2_const")
  Case_2_const = fit_and_save("Case_2_const", Case_2_stan_data)
  Case_2_const_bridge = bridgesampling::bridge_sampler(Case_2_const, silent = TRUE)
  saveRDS(Case_2_const_bridge,
          file = file.path("Model Fitting", "Fit Models","Case_2_const_bridge.rds"))
  
  log_progress("[6/9] Fitting Case_2_true")
  Case_2_true = fit_and_save("Case_2_true", Case_2_stan_data)
  Case_2_true_bridge = bridgesampling::bridge_sampler(Case_2_true, silent = TRUE)
  saveRDS(Case_2_true_bridge,
          file = file.path("Model Fitting", "Fit Models","Case_2_true_bridge.rds"))
  
} # Case 2 

suppressWarnings({ # The jumping between positive and negative biases creates a problematic fit, hence the warning suppression
  log_progress("[7/9] Fitting Case_2_extra_random")
  Case_2_extra_random         = fit_and_save("Case_2_extra_asymm", Case_2_extra_random_data,   suffix = "_random")
  
  log_progress("[8/9] Fitting Case_2_extra_opposing")
  Case_2_extra_opposing       = fit_and_save("Case_2_extra_asymm", Case_2_extra_opposing_data, suffix = "_opposing")
  
  log_progress("[9/9] Fitting Case_2_extra_opposing_non_H")
  Case_2_extra_opposing_non_H = fit_and_save("Case_2_extra_asymm_non_H", Case_2_extra_opposing_data, suffix = "_opposing")
}) # Case 2 extras
