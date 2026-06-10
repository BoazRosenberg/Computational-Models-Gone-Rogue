# Computational Models Gone Rogue

Code accompanying the manuscript *Computational Models Gone Rogue*.

## Overview

This repository contains the code used to generate the simulated datasets, fit the computational models, perform posterior predictive checks and additional simulations, and reproduce the analyses reported in the manuscript.

## Reproducing the Results

### Step 1: Generate the Simulated Data

Run:

`simulate_data.py`

located in the `Data Simulation` directory.

This script generates all datasets used in the manuscript and saves them to:

`Data Simulation/Simulated Data`

### Step 2: Run the Analysis Pipeline

Open `Gone Rogue.Rproj` and run:

`source("run_all.R")`

The pipeline will:

1. Fit all Stan models.
2. Run posterior predictive checks and additional simulations.
3. Render the results report.

## Repository Structure

### Root Directory

Key files for reproducing the analyses:

* `Gone Rogue.Rproj` is the RStudio project file for the repository.
* `run_all.R` executes the full analysis pipeline, including model fitting, simulations, and report generation.
* `report_results.Rmd` renders the final analysis report presented in the manuscript.

### Data Simulation

Python code used to generate the simulated datasets analyzed in the manuscript.

* `functions.py` contains the task-data-simulation function `create_data`.
* `models.py` contains the generative models used to simulate choice data.
* `simulate_data.py` generates all datasets for Cases 1 and 2, and extra variants of Case 2 with a two-armed bandit task with full information.

### Model Fitting

R and Stan code used to fit the computational models to the simulated datasets.

* `fit_stan_models.R` fits all models.
* `setup.R` contains supporting functions and setup code.
* `Stan Models/` contains the Stan model specifications.
* `Fit Models/` contains saved model fits and bridge-sampling results.

### Simulations and Analysis

Scripts used to generate the analyses reported in the manuscript.

* `PPC1.R` performs posterior predictive checks for Case 1.
* `PPC2.R` performs posterior predictive checks for Case 2.
* `simulate_Case_2.R` simulates model update dynamics for Case 2.
* `RDS Files/` contains saved simulation outputs.

### Figures

All figures produced in `report_results.Rmd` are automatically exported as PNG files during rendering.

These files correspond to the figures reported in the manuscript.

## Notes

* The Python simulation step must be completed before running the R pipeline.
* Individual scripts can be executed separately, but `run_all.R` provides the recommended workflow for reproducing all analyses and results.
* Model fitting is computationally intensive and may require substantial runtime and memory resources.

## Citation

If you use this code, please cite:

Rosenberg, B. & Eldar E. (2026). *Computational Models Gone Rogue*.