# ==============================================================================
# Epidemic simulation QC script
# ==============================================================================
#
# PURPOSE
# -------
#
# APPROACH
# --------
# Loads output from simulation.R
# Benchmarking guidance mainly from:
# https://academic.oup.com/ve/article/9/1/vead010/7028398#401056156
# and https://onlinelibrary.wiley.com/doi/10.1002/sim.8086
#
# OUTPUTS
# -------
#
# ============================================================================== 


# Validation should occur on two leves:
# Within each simulation, and across all simulations

# Tests to do:
# - Root-to-tree regression
# - Compare R0

library(ggplot2)
library(dplyr)


simulation_dirs <- list.dirs("output/simulation/")

# Remove first and last elements
simulation_dirs <- simulation_dirs[-1]  # Root folder
simulation_dirs <- simulation_dirs[-length(simulation_dirs)] # log folder

# For now, work with only one, but later create a loop
path <- simulation_dirs[1]

sim <- readRDS(paste0(path, "/nosoi_sim.rds"))

host_table <- getTableHosts(sim)
state_table <- getTableState(sim)
dynamics_table <- getDynamic(sim)
cumu_table <- getCumulative(sim)
sim_sum <- summary(sim)


ggplot(cumu_table, aes(t, Count)) +
  geom_line() 

ggplot(dynamics_table |> filter(state == "Netherlands"), aes(t, Count, color = state)) +
  geom_line()


# Load all dynamics tables into one dataframe
i <- as.integer(1)
dynamics_all <- data.frame()

for (dir in simulation_dirs) {
  
  file <- paste0(dir, "/nosoi_sim.rds")
  
  if (file.exists(file)) {
    sim <- readRDS(file)
    
    dynamics_table <- getDynamic(sim)
    dynamics_table$sim_id <- toString(i)
    
    dynamics_all <- rbind(dynamics_all, dynamics_table)  
  } else {
    message(paste0("Simulation missing in ", file))
  }
  
  i <- i + 1L
}

dynamics_avg <- dynamics_all |> 
  group_by(t, state) |> 
  summarise(Count = mean(Count, na.rm = TRUE), .groups = "drop")

ggplot() +
  # Individual simulations: faded, grouped so lines don't connect across sim_id
  geom_line(data = dynamics_all, 
            aes(t, Count, group = sim_id), 
            color = "grey70", alpha = 0.4, linewidth = 0.3) +
  # Average trajectory on top
  geom_line(data = dynamics_avg, 
            aes(t, Count), 
            color = "firebrick", linewidth = 1) +
  facet_wrap(~state) +
  theme_minimal()
