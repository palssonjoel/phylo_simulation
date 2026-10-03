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


simulation_dirs <- list.dirs("output/simulation_test_density/")

# Remove first and last elements
simulation_dirs <- simulation_dirs[-1]  # Root folder
simulation_dirs <- simulation_dirs[-length(simulation_dirs)] # log folder

############################3
# TESTING
path <- simulation_dirs[1]

sim <- readRDS(paste0(path, "1/nosoi_sim.rds"))

host_table <- getTableHosts(sim)
state_table <- getTableState(sim)
dynamics_table <- getDynamic(sim)
cumu_table <- getCumulative(sim)
sim_sum <- summary(sim)

ggplot(cumu_table, aes(t, Count)) +
  geom_line() 

ggplot(dynamics_table |> filter(state == "Denmark"), aes(t, Count, color = state)) +
  geom_line() +
  labs(
    title = "No. active hosts"
  )

# R0
data = data.frame(R0=getR0(sim)$R0.dist)
ggplot(data=data, aes(x=R0)) + geom_histogram() + theme_minimal()

# Active hosts per timestep, by state
dynamics_table |> 
  group_by(state, t) |> 
  ggplot(aes(t, Count, color=state)) +
  geom_line() +
  labs(
    title = "Epidemic dynamics",
    y = "No. Active hosts",
    x = "time"
  )

# no introductions into states
# Look up the infection state of each infector
infector_states <- host_table[, .(
  infector_id = hosts.ID,
  infector_state = inf.in
)]

# Add infector state to each infected host
host_table2 <- merge(
  host_table,
  infector_states,
  by.x = "inf.by",
  by.y = "infector_id",
  all.x = TRUE
)

introductions <- host_table2[
  !is.na(infector_state) &
    infector_state != inf.in,
  .(introductions = .N),
  by = inf.in
]

setorder(introductions, -introductions)

introductions

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
  geom_line(data = test$dynamics_all, 
            aes(t, Count, group = sim_id), 
            color = "grey70", alpha = 0.4, linewidth = 0.3) +
  # Average trajectory on top
  geom_line(data = dynamics_avg, 
            aes(t, Count), 
            color = "firebrick", linewidth = 1) +
  facet_wrap(~state) +
  ylim(0, 100) +
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



##########################
#TESTING BELOW
# Animate a nosoi discrete-structure simulation across European countries
# Adapted from https://slequime.github.io/nosoi/articles/examples/viz.html
# (the "Benelux" map example), rewritten to scale to 18 states.
#
# Install once if needed:
# install.packages(c("data.table","dplyr","ggplot2","gganimate","viridis","maps","gifski"))

library(nosoi)
library(data.table)
library(dplyr)
library(ggplot2)
library(gganimate)
library(viridis)

# ---- 0. Settings -------------------------------------------------------------
sim    <- sim        # <-- replace with the name of YOUR nosoiSim object
metric <- "active"           # "active"     = hosts currently infected in each country
# "cumulative" = total infections that have occurred there
step   <- 1L                 # time-step thinning: 1 = every step, 5 = every 5th step, etc.
# (raise this if your simulation is long; frames = total.time/step)

# ---- 1. Country centroids (state names must match your simulation's states) ---
centroids <- data.table(
  name = c("Albania","Austria","Belgium","Croatia","Czechia","Denmark","Germany",
           "Greece","Hungary","Italy","Luxembourg","Netherlands","Poland",
           "Portugal","Romania","Serbia","Slovakia","Spain"),
  lat  = c(41.15, 47.52, 50.50, 45.10, 49.82, 56.26, 51.17,
           39.07, 47.16, 41.87, 49.82, 52.13, 51.92,
           39.40, 45.94, 44.02, 48.67, 40.46),
  long = c(20.17, 14.55,  4.47, 15.20, 15.47,  9.50, 10.45,
           21.82, 19.50, 12.57,  6.13,  5.29, 19.15,
           -8.22, 24.97, 21.01, 19.70, -3.75)
)

centroids <- rbind(centroids, data.table(
  name = c("Bulgaria", "France", "Slovenia"),
  lat  = c(42.73, 46.23, 46.15),
  long = c(25.49,  2.21, 14.99)
))

# ---- 2. Extract state occupancy over time -------------------------------------
st   <- as.data.table(getTableState(sim))   # hosts.ID, state, time.from, time.to
tmax <- as.integer(sim$total.time)

st[, `:=`(time.from = as.integer(time.from), time.to = as.integer(time.to))]
st[is.na(time.to), time.to := tmax + 1L]    # still-active hosts stay until the end
st[, row := .I]

frames <- seq(0L, tmax, by = step)

# Which state each host is in at each time step (vectorised replacement for the
# slow for-loop used in the tutorial)
occupancy <- st[, .(time = time.from:max(time.from, time.to - 1L)), by = row]
occupancy[, state := st$state[row]]
occupancy <- occupancy[time %in% frames]

active_counts <- occupancy[, .N, by = .(time, state)]

# Cumulative infections per country, from the infection table (inf.in / inf.time)
# host_table must be your getTableHosts() output, as in your message
hosts <- as.data.table(host_table)[!is.na(inf.in)]
hosts[, inf.time := as.integer(inf.time)]
cum_counts <- CJ(state = centroids$name, time = frames)
cum_counts[, N := mapply(function(s, t) sum(hosts$inf.in == s & hosts$inf.time <= t),
                         state, time)]

counts <- if (metric == "active") active_counts else cum_counts
counts <- counts[N > 0] %>%
  left_join(centroids, by = c("state" = "name"))

# ---- 3. Movement events between countries ------------------------------------
setorder(st, hosts.ID, time.from)
moves <- st[, .(from = state, to = shift(state, -1L), time = time.to), by = hosts.ID
][!is.na(to) & from != to]
moves[, time := (time %/% step) * step]     # snap to the animation frames
moves <- moves[time %in% frames, .N, by = .(time, from, to)]

moves <- moves %>%
  left_join(centroids, by = c("from" = "name")) %>%
  rename(lat_from = lat, long_from = long) %>%
  left_join(centroids, by = c("to" = "name")) %>%
  rename(lat_to = lat, long_to = long)

# Catch name mismatches early
stopifnot("Some states in the simulation don't match `centroids$name`" =
            !anyNA(counts$lat) && !anyNA(moves$lat_from) && !anyNA(moves$lat_to))

# ---- 4. Animated map ---------------------------------------------------------
worldmap <- borders("world", colour = "gray50", fill = "#efede1")

anim_plot <- ggplot() +
  worldmap +
  coord_fixed(ratio = 1.5, xlim = c(-10, 30), ylim = c(35, 58)) +
  # movement arrows (thickness = number of movers in that step)
  geom_curve(data = moves,
             aes(x = long_from, y = lat_from, xend = long_to, yend = lat_to,
                 linewidth = N),
             arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
             curvature = 0.2, colour = "gray30", alpha = 0.6) +
  # bubbles per country
  geom_point(data = counts,
             aes(long, lat, colour = state, size = N), alpha = 0.8) +
  # country labels (static)
  geom_text(data = centroids, aes(long, lat - 0.9, label = name),
            size = 2.5, colour = "black") +
  scale_size_area(max_size = 25, guide = "none") +
  scale_linewidth(range = c(0.3, 1.5), guide = "none") +
  scale_colour_viridis_d(guide = "none") +
  theme_void() +
  transition_states(time, transition_length = 1, state_length = 1) +
  labs(title = paste0("Time: {closest_state}  (", metric, " infections)"))

anim <- animate(anim_plot,
                nframes = length(frames) * 2 + 10,
                fps = 10, width = 900, height = 800, end_pause = 10)

anim_save("nosoi_europe.gif", anim)

# For an mp4 instead (needs the 'av' package):
# animate(anim_plot, nframes = length(frames) * 2 + 10, fps = 10,
#         width = 900, height = 800, renderer = av_renderer("nosoi_europe.mp4"))