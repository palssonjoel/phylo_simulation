# ==============================================================================
# Subsampling scheme
# ==============================================================================
#
# PURPOSE
# -------
#
# APPROACH
# --------
#
#
# OUTPUTS
# -------
#
# ============================================================================== 

# Three proposed schemes:
# 1) Even sampling across all demes
# 2) Sampling proportional to deme population
# 3) Biased sampling mimicing real life

# Even sampling should be easy, question is: is it an even proportion (e.g. 
# 20% of all cases in every deme), or an even number (e.g. 50 samples in each deme)?

# Proportional is also simple. Calculate population proportion, and apply to
# the entire set of hosts to be sampled. If Denmark has 50% of animals,
# they have 50% of cases in the subsampled data.

# Biased sampling is built on same approach, but proportions are instead calculated
# based on proportion of sampling in datasets.

library(dplyr)

source("code/simulation.R")

simulation_data <- read.csv("output/simulation/1/simulation_data.csv")
pop_stats <- read.csv("output/trade_movement_rates_eu.csv") |> 
  select(exporter, exporter_herd_size) |> 
  distinct()

# Calculate proportion of population
sampling_proportions          <- pop_stats
sampling_proportions$total    <- sum(sampling_proportions$exporter_herd_size)
sampling_proportions$pop_prop <- sampling_proportions$exporter_herd_size / sampling_proportions$total 
sampling_proportions |> select(exporter, pop_prop)

# For CPU hour estimate, a random sampling scheme will be applied

# Approach:
# 1. Load the nosoi + hky output
# 2. Extract rows in a random manner
# 3. Apply create_alignment() from simulation.R
# 4. Done! :) 

random_subsample <- slice_sample(simulation_data, n=150)
create_alignment(random_subsample, "random_subsample", "output/simulation/1")
