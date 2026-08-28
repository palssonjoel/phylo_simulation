# ==============================================================================
# Main swIAV epidemic simulatuion
# ==============================================================================
#
# PURPOSE
# -------
#
# APPROACH
# --------
# Based on the tutorial: https://slequime.github.io/nosoi/articles/discrete.html
#
# OUTPUTS
# -------
#
# ============================================================================== 

library(nosoi)
library(tidyverse)

# Load population data for max infections
pop_data <- read_csv("output/trade_movement_rates.csv") |> 
  select(exporter, exporter_herd_size) |> 
  distinct(exporter, exporter_herd_size) |> 
  rename(country = exporter,
         swine_heads = exporter_herd_size)

# Load transition matrix
transition_matrix <- as.matrix(
  read.csv("output/trade_matrix_daily_probabilities.csv",
           row.names = 1,
           check.names = FALSE)
)

# Normalize to avoid floating point issue (twice on purpose)
transition_matrix <- transition_matrix / rowSums(transition_matrix)
transition_matrix <- transition_matrix / rowSums(transition_matrix)

rowSums(transition_matrix) - 1
# Set simulation functions

# Core  set of functions for the simulation are:
#   - pExit = Probability that host exits simulation (death, cured, etc.)
#   - pMove = Probability that a host leaves its current state. NOT the same 
#             as movement probabilities in transition matrix.
#   - sdMove = SD of the random walk in pMove
#   - nContact = Number of potentially infectious contacts an infected host can 
#                encounter per unit of time.
#   - pTrans = Probability of transmission over time, when a contact occurs. Can
#              be set to e.g. seasionality. 
# 

# pExit is assumed to equal death rate due to swIAV. As per one article, it is
# between 10-15%, with up to 100% morbidity:
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7587018/

p_exit_func <- function(t) {
  runif(1, min = 0.10, max = 0.15) # random value between 10-15%
}

# pMove is probability that one pig moves between locations. Settings to 10% for now
p_move_func <- function(t) { 
  return(0.1)
}

# nContact is not equal to R0 (number of 2nd infections) - pContact is only how
# many individuals a host encounters. This is HIGHLY dependent on farm structure
# age of pig, whether it's transported etc. Need to find data on this.
# I will infer it from pig pen sizes. According to article below, it is 12 pigs 
# per pen in Europe. There is likely variation to this
# https://ahdb.org.uk/eupig-novelty-in-enrichment-material

# ASSUMPTION: For now, I will assume that it is the same for all countries.
n_contact_func <- function(t) {
  abs(round(rnorm(1, 12, 3)))
}

# pTrans depends on incubation time, which is often cited between 1-3 days
# p_max is a constant probability of transmission
p_trans_func <- function(t, p_max, t_incubation) {
  if(t < t_incubation){p = 0}
  if(t > t_incubation){p = p_max}
  return(p)
}

t_incub_func <- function(x){pmax(0, rnorm(x, mean = 2, sd = 0.5))}
p_max_func <- function(x){rbeta(x, shape1=5, shape2=2)}

################################################################################
# Epidemic dynamics histograms
################################################################################

# Number of simulated values
n <- 10000

# Generate values
t_incub <- t_incub_func(n)
p_max <- p_max_func(n)
n_contacts <- replicate(n, n_contact_func(0))

# Histograms
hist(
  t_incub,
  breaks = 50,
  main = "Distribution of incubation time",
  xlab = "Incubation time",
  ylab = "Frequency"
)

hist(
  p_max,
  breaks = 50,
  main = "Distribution of maximum transmission probability",
  xlab = "p_max",
  ylab = "Frequency"
)

hist(
  n_contacts,
  breaks = seq(-0.5, max(n_contacts) + 0.5, by = 1),
  main = "Distribution of nContact",
  xlab = "Number of contacts",
  ylab = "Frequency"
)

################################################################################
# Run simulation
################################################################################
set.seed(1998)
max_infections = 10000

SimulationSingle <- nosoiSim(type="single", popStructure="discrete",
                             length.sim=1095, max.infected=max_infections, init.individuals=1, init.structure="Denmark", 
                             
                             structure.matrix=transition_matrix,
                             
                             pExit = p_exit_func,
                             param.pExit=NA,
                             timeDep.pExit=FALSE,
                             diff.pExit=FALSE,
                             
                             pMove = p_move_func,
                             param.pMove=NA,
                             timeDep.pMove=FALSE,
                             diff.pMove=FALSE,
                             
                             nContact=n_contact_func,
                             param.nContact=NA,
                             timeDep.nContact=FALSE,
                             diff.nContact=FALSE,
                             
                             pTrans = p_trans_func,
                             param.pTrans = list(p_max=p_max_func,t_incubation=t_incub_func),
                             timeDep.pTrans=FALSE,
                             diff.pTrans=FALSE,
                             
                             prefix.host="H",
                             print.progress=TRUE,
                             print.step=10)











