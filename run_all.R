# =========================================================
# Master Pipeline Script
# =========================================================
#
# IMPORTANT (PYTHON PREPROCESSING):
# Before running this R pipeline, the Python simulation step must be executed.
#
# The Python code is located in:
#   Data Simulation/
#
# It contains three scripts:
#   - functions.py
#   - models.py
#   - simulate_data.py
#
# Only the following script needs to be executed:
#   simulate_data.py
#
# This script generates the simulated data sets we will be fitting our Stan model to

# ----------------------------------------
# Helper function for reporting progress
# ----------------------------------------

{
start_time <- Sys.time()

log_progress <- function(msg) {
  now <- Sys.time()
  elapsed <- as.numeric(difftime(now, start_time, units = "mins"))
  cat(sprintf("[%s | +%.1f min] %s\n", format(now, "%H:%M"), elapsed, msg),
      file = stderr())
}

run_step <- function(script_path) {

  # Create a temporary sandbox environment and run script
  temp_env <- new.env()
  temp_env$start_time <- start_time
  
 
    tryCatch({
      source(script_path, local = temp_env)
    }, error = function(e) {
      log_progress(sprintf("! Error in %s: %s\n", script_path, e$message))
    })
 
  
  # Cleanup: Delete the sandbox and force garbage collection
  rm(temp_env)
  invisible(gc(verbose = FALSE))

}
}

# ---------------------------
# 1. Fit all Stan models  
# ---------------------------

# First, we fit all Stan models to the corresponding data sets.

# This may take a substantial amount of computing power and memory, 
# as it includes fitting 9 Stan models sequentially.

# Each model is fit using:
#   - 4 MCMC chains
#   - 2,500 iterations per chain
{
log_progress("Starting Step 1: Model Fitting")
run_step("Model Fitting/fit_stan_models.R")
}
# ---------------------------
# 2. Run Simulations and Posterior Predictive Checks
# ---------------------------

# These scripts generate outputs that are saved to disk (RDS files).
# They are less computationally intensive than model fitting,
# but may still require substantial RAM.
{
cat("\n")
log_progress("Starting Step 2: Simulations and Posterior Predictive Checks")
log_progress("[1/3] Simulating PPC for Case 1")
run_step("Simulations and Analysis/PPC1.R")

log_progress("[2/3] Simulating PPC for Case 2")
run_step("Simulations and Analysis/PPC2.R")

log_progress("[3/3] Simulating mean updates for Case 2")
run_step("Simulations and Analysis/simulate_Case_2.R")
}

# ---------------------------
# 3. Knit R Markdown report
# ---------------------------

# Knit report
{
cat("\n")
log_progress("Rendering R Markdown report")
out = rmarkdown::render(
  input = "report_results.Rmd",
  quiet = T
)

log_progress("Pipeline completed successfully.")

log_progress(sprintf("Results report created: %s", out))
browseURL(out)
}

