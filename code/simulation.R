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
library(seqinr)
library(phangorn)

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

# pMove is probability that one pig moves between locations. Settings to 10% for now
p_move_func <- function(t) { 
  return(0.1)
}

# nContact is not equal to R0 (number of 2nd infections) - pContact is only how
# many individuals a host encounters. This is HIGHLY dependent on farm structure
# age of pig, whether it's transported etc. Need to find data on this.
# I will infer it from pig pen sizes. According to article below, it is 12 pigs 
# per pen in Europe. There is likely variation to this. 
# However, in a pen, 12 pigs doesn't mean 12 naive hosts, so nContacts can't equal this.
# If it does, it causes explosive growth. I've adjusted it down. 
# https://ahdb.org.uk/eupig-novelty-in-enrichment-material

# ASSUMPTION: For now, I will assume that it is the same for all countries.
n_contact_func <- function(t) {
  abs(round(rnorm(1, 2, 1)))
}

# pTrans depends on incubation time, which is often cited between 1-3 days
# p_max is a constant probability of transmission
p_trans_func <- function(t, p_max, t_incubation) {
  if(t < t_incubation){p = 0}
  if(t > t_incubation){p = p_max}
  return(p)
}

t_incub_func <- function(x){pmax(0, rnorm(x, mean = 2, sd = 0.5))}
p_max_func <- function(x){rbeta(x, shape1=1, shape2=3)}


# As per one article,  death rate due to swIAV is
# between 10-15%, with up to 100% morbidity. But this is over the course of an
# entire illness, I'll use a duration based exit, so hosts can't exit before
# the incubation period is done.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7587018/

# Approach: zero exit probability duing incubation, then a 1/5 probability of exiting
# per day. 

p_exit_func <- function(t, t_incubation) {
  if (t < t_incubation) { return(0) }
  else {
    return(1/5)  # ≈ 0.20/day → mean ~5 days post-incubation illness
  }
}

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
# Run transmission simulation
################################################################################
set.seed(12092)
max_infections = 10000

simulation <- nosoiSim(type="single", popStructure="discrete",
                             length.sim=365, max.infected=max_infections, init.individuals=1, init.structure="Denmark", 
                             
                             structure.matrix=transition_matrix,
                             
                             pExit = p_exit_func,
                             param.pExit=list(t_incubation=t_incub_func),
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


tree <- getTransmissionTree(simulation)


ggtree(tree) + 
  geom_nodepoint(aes(color=state)) + 
  geom_tippoint(aes(color=state)) +
  theme_tree2() + xlab("Time (t)") + theme(legend.position = c(0.05,0.8), 
                                           legend.title = element_blank(),
                                           legend.key = element_blank())



################################################################################
# HKY model
################################################################################
# Conceptually, to model molecular evolution, I think I only need:
# - Who infected who
# - when they were infected
# This gives  duration of infection (in days from nosoi), to which I can apply 
# a substitution rate. I have all this informatoin from getTableHosts():
# - host.ID = ID of the host
# - inf.by = who the host was infected by
# - inf.time = the timepoint host was infected (entered simulation)
# - out.time = the timepoint the host left the simulation

calculate_inf_duration <- function(host_table) {
  # Calculate the evolutionary time to use for the HKY substition model.
  # Evolutionary time is calculated as:
  # infector's time of infect - infectee's time of infection
  
  # For each host, find the row index of their infector in host_table
  infector_idx <- match(host_table$inf.by, host_table$hosts.ID)
  
  # Look up the infector's inf.time using those indices
  infector_time <- host_table$inf.time[infector_idx]
  
  # Evolution time = own inf.time minus infector's inf.time
  host_table$evo.time <- host_table$inf.time - infector_time 
  
  return(host_table)
}

# HKY ALGORITHM
hky_nosoi <- function(ref_genome, mu = 7.6e-3, kappa = 4.5, host_data) {
  # The algorithm:
  # For each host (except index case), extract the sequence of who infected them.
  # Then, apply a HKY substitution model to this sequence.
  # HKY is a function of time, meaning that substitution rates occur with a probability
  # over time. Here, what is relevant is the difference between when the thost was infected
  # and when their infector was infected. 
  
  # Variables: 
  # ref_genome = reference genome, that of the virus infecting first host
  # mu = yearly substition rate, default 7.6e-3
  # kappa = transition/transversion rate ratio, default 4.5
  # host_data = nosoi getHostData() output
  
  # HKY PARAMETERS
  bases <- c("a", "c", "g", "t")
  mu <- mu / 365.25 # HA mutation rate per day
  baseComp <- table(ref_genome)
  baseFreq <- as.numeric(baseComp) / length(ref_genome)
  names(baseFreq) <- names(baseComp)
  
  # beta is "base rate" appied to transversions
  # Formula explained: given these base frequencies and this kappa, 
  # what value of beta makes the overall average substitution rate come out to exactly mu?
  beta <- mu/(2*(baseFreq['a']*baseFreq['c']+baseFreq['a']*baseFreq['t']+baseFreq['c']*baseFreq['g']+baseFreq['g']*baseFreq['t'])+
                2*kappa*(baseFreq['a']*baseFreq['g']+baseFreq['c']*baseFreq['t']))
  
  # alpha is transtion rate
  alpha <- kappa * beta
  
  # Q matrix
  # In HKY, mutations are dependent on whether it's a transition or transversion
  # and how common the destination base is, hence *baseFre
  Q <- matrix(NA, ncol = 4, nrow = 4)
  Q[upper.tri(Q)] <- c(beta*baseFreq['c'], alpha*baseFreq['g'], 
                       beta*baseFreq['g'], beta*baseFreq['t'], 
                       alpha*baseFreq['t'], beta*baseFreq['t']) 
  Q[lower.tri(Q)] <- c(beta*baseFreq['a'], alpha*baseFreq['a'], 
                       beta*baseFreq['a'], beta*baseFreq['c'], 
                       alpha*baseFreq['c'], beta*baseFreq['g'])
  diag(Q) <- -apply(Q, 1, sum, na.rm = TRUE)
  
  # Compute the invertible matrix of Q and its eigenvalues 
  # This is done only once to save compuational load
  eig <- eigen(Q)
  E <- eig$vectors       # matrix of eigenvectors
  eigvals <- eig$values  # vector of eigenvalues
  E_1 <- solve(E)        # inverse of E
  
  # Construct host_table
  host_table <- host_data |> 
    select(hosts.ID, inf.by, inf.time, out.time) 
  
  host_table <- calculate_inf_duration(host_table) # Calc. evolutionary time
  
  # Add index sequence
  host_table$seq <- vector("list", nrow(host_table))
  host_table$seq[[1]] <- ref_genome # Apend ref genome to first host
  
  # Apply HKY to sequences
  evolve <- function(host_table, E, E_1, eigvals, bases) {
    
    # Create a list of sequences for all hosts
    n <- nrow(host_table)
    seq_list <- vector("list", n)
    names(seq_list) <- host_table$hosts.ID
    seq_list["H-1"] <- host_table$seq[1]  # First host get reference genome
    
    host_nr <- order(host_table$inf.time) # Extract host # ordered by inf.time
    host_nr <- setdiff(host_nr, 1)        # Skip first host, already done
    
    # Process all hosts in order of infection
    for (i in host_nr) {
      infector_id <- host_table$inf.by[i]
      t <- host_table$evo.time[i]
      
      # %*% = matrix multiplication operator
      P_t <- Re(E %*% diag(exp(eigvals * t)) %*% E_1)
      dimnames(P_t) <- list(bases, bases)
      
      infector_seq <- seq_list[[infector_id]]
      
      # Apply substitution independently at each site
      new_seq <- vapply(infector_seq,
                        function(b) {
                          # Sample a base with the probability of P_t
                          sample(bases, size = 1, prob = P_t[b, ])
                        }, character(1))
      
      seq_list[[host_table$hosts.ID[i]]] <- new_seq
    }
    
    return(seq_list)
  }
  
  sequences <- evolve(host_table, E, E_1, eigvals, bases)
  host_table$seq <- sequences[host_table$hosts.ID]
  
  return(host_table)
  
}


host_data <-  getHostData(simulation)
ref_genome <- read.fasta("data/ref_genome.fasta", forceDNAtolower = FALSE, set.attributes = FALSE)[[1]]

test <- hky_nosoi(ref_genome = ref_genome, host_data = host_data)
