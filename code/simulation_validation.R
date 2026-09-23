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
library(nosoi)


simulation_dirs <- list.dirs("output/simulation/")

# Remove first and last elements
simulation_dirs <- simulation_dirs[-1]  # Root folder
simulation_dirs <- simulation_dirs[-length(simulation_dirs)] # log folder

############################3
# TESTING
path <- simulation_dirs[1]

sim <- readRDS(paste0(path, "/nosoi_sim.rds"))

host_table <- getTableHosts(sim)
state_table <- getTableState(sim)
dynamics_table <- getDynamic(sim)
cumu_table <- getCumulative(sim)
sim_sum <- summary(sim)


ggplot(cumu_table, aes(t, Count)) +
  geom_line() 

ggplot(dynamics_table |> filter(state == "Denmark"), aes(t, Count, color = state)) +
  geom_line()

#########################################################
# Load all dynamics tables into one dataframe

generate_state_graphs <- function(simulation_dirs) {
  
  host_all     <- data.frame()
  state_all    <- data.frame()
  dynamics_all <- data.frame()
  cumu_all     <- data.frame()
  
  i <- as.integer(1) 
  
  for (dir in simulation_dirs) {
    sim_file <- paste0(dir, "/nosoi_sim.rds")
    
    if (file.exists(sim_file)) {
      sim <- readRDS(sim_file)
      
      # Load data tables
      host_table     <- getTableHosts(sim)
      state_table    <- getTableState(sim)
      dynamics_table <- getDynamic(sim)
      cumu_table     <- getCumulative(sim)
      
      # Add simulation ID
      host_table$sim_id     <- toString(i)
      state_table$sim_id    <- toString(i)
      dynamics_table$sim_id <- toString(i)
      cumu_table$sim_id     <- toString(i)
      
      
      host_all     <- rbind(host_all, host_table)
      state_all    <- rbind(state_all, state_table)  
      dynamics_all <- rbind(dynamics_all, dynamics_table)  
      cumu_all     <- rbind(cumu_all, cumu_table)  
    } else {
      message(paste0("Simulation missing in ", sim_file))
    }
    i <- i + 1L
  }
  
  all_data <- list(
    host_all     = host_all,
    state_all    = state_all,
    dynamics_all = dynamics_all,
    cumu_all     = cumu_all
  )
  
  return(all_data)
}

test <- generate_state_graphs(simulation_dirs)

dynamics_avg <- test$dynamics_all |> 
  group_by(t, state) |> 
  summarise(Count = mean(Count, na.rm = TRUE),
            sd = sd(Count, na.rm = TRUE),
            .groups = "drop")

cumu_avg <- test$cumu_all |> 
  group_by(t) |> 
  summarise(Count = mean(Count, na.rm = TRUE),
            sd = sd(Count, na.rm = TRUE),
            .groups = "drop")

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

# Cumulative 
ggplot() +
  # Individual simulations: faded, grouped so lines don't connect across sim_id
  geom_line(data = test$cumu_all, 
            aes(t, Count, group = sim_id), 
            color = "grey70", alpha = 0.4, linewidth = 0.3) +
  # Average trajectory on top
  geom_line(data = cumu_avg, 
            aes(t, Count), 
            color = "firebrick", linewidth = 1) +
  theme_minimal()
