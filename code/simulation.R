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
library(dplyr)
library(ggtree)
library(ggplot2)
library(treeio)
library(seqinr)
library(phangorn)
library(logr)
library(Biostrings)

run_nosoi <- function(transition_matrix, max_infections, simulation_time, out_dir, sim_length) {

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
 incub_hist <- hist(
    t_incub,
    breaks = 50,
    main = "Distribution of incubation time",
    xlab = "Incubation time",
    ylab = "Frequency"
  )
  
 p_max_hist <- hist(
    p_max,
    breaks = 50,
    main = "Distribution of maximum transmission probability",
    xlab = "p_max",
    ylab = "Frequency"
  )
  
  n_contacts_hist <- hist(
    n_contacts,
    breaks = seq(-0.5, max(n_contacts) + 0.5, by = 1),
    main = "Distribution of nContact",
    xlab = "Number of contacts",
    ylab = "Frequency"
  )
  
  # Save plots
  png(file.path(out_dir, "incub_hist.png"), width = 800, height = 600)
  plot(incub_hist, main = "Distribution of incubation time", xlab = "Incubation time")
  dev.off()
  
  png(file.path(out_dir, "p_max_hist.png"), width = 800, height = 600)
  plot(p_max_hist, main = "Distribution of maximum transmission probability", xlab = "p_max")
  dev.off()
  
  png(file.path(out_dir, "n_contacts_hist.png"), width = 800, height = 600)
  plot(n_contacts_hist, main = "Distribution of nContact", xlab = "Number of contacts")
  dev.off()
  
  ################################################################################
  # Run transmission simulation
  ################################################################################
  
  msg <- (paste0("Starting transmission simulation\n",
                   "Timepoint: ", format(Sys.time(), "%H:%M:%S"), "\n",
                   "Max infected: ", max_infections, "\n",
                   "Time limit: ", sim_length, "\n"))
  log_print(msg)

  # Save nosoi simulation output in log
  lf <- log_path()
  con <- file(lf, open = "a")
  sink(con, split = TRUE, type = "output")
  
  # Run simulation
  time_start <- Sys.time()
  simulation <- nosoiSim(type="single", popStructure="discrete",
                         length.sim=sim_length, max.infected=max_infections, init.individuals=1, init.structure="Denmark", 
                         
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
  time_end <- Sys.time()
  duration <- difftime(time_end, time_start, units = "mins")
  
  sink(type = "output")
  close(con)
  
  msg <- paste0(
    "Transmission simulation complete. Time elapsed: ",
    round(as.numeric(duration), 2), " ", units(duration)
  )
  log_print(msg)
  
  log_print("=========================================================================")
  
  # Get and save transmission tree
  tree <- getTransmissionTree(simulation)
  write.beast(tree, file.path(out_dir, "transmission_tree.nexus"))
  write.beast.newick(tree, file.path(out_dir, "transmission_tree.nwk"))
 
  return(simulation)
}

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
  mu_daily <- mu / 365.25 # HA mutation rate per day
  baseComp <- table(ref_genome)
  baseFreq <- as.numeric(baseComp) / length(ref_genome)
  names(baseFreq) <- names(baseComp)
  
  # beta is "base rate" appied to transversions
  # Formula explained: given these base frequencies and this kappa, 
  # what value of beta makes the overall average substitution rate come out to exactly mu?
  beta <- mu_daily/(2*(baseFreq['a']*baseFreq['c']+baseFreq['a']*baseFreq['t']+baseFreq['c']*baseFreq['g']+baseFreq['g']*baseFreq['t'])+
                2*kappa*(baseFreq['a']*baseFreq['g']+baseFreq['c']*baseFreq['t']))
  
  # alpha is transion rate
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
 # host_table <- host_data |> 
  #  select(hosts.ID, inf.by, inf.time, out.time) 
  
  host_table <- calculate_inf_duration(host_data) # Calc. evolutionary time
  
  # Add index sequence
  host_table$seq <- vector("list", nrow(host_table))
  host_table$seq[[1]] <- ref_genome # Apend ref genome to first host
  
  # Make readable strings for log
  base_comp_str <- paste(names(baseComp), baseComp, sep = "=", collapse = ", ")
  base_freq_str <- paste(names(baseFreq), round(baseFreq, 4), sep = "=", collapse = ", ")
  
  msg <- paste0(
    "Applying HKY substitution model\n",
    "Timepoint: ", format(Sys.time(), "%H:%M:%S"), "\n\n",
    "Mutation rate: ", mu, "\n",
    "Transversion rate: ", beta, "\n",
    "Transition rate: ", alpha, "\n",
    "Kappa: ", kappa, "\n\n",
    "Reference base composition: ", base_comp_str, "\n")
  log_print(msg)
  
  time_start <- Sys.time()
  
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
  
  time_end <- Sys.time()
  duration <- difftime(time_end, time_start, units = "mins")
  
  msg <- paste0(
    "HKY substitution complete. Time elapsed: ",
    round(as.numeric(duration), 2), " ", units(duration)
  )
  
  log_print(msg)
  
  return(host_table)
}

run_simulation <- function(max_infections, 
                           sim_length,
                           seed, 
                           out_dir) {
  # This runs both the transmission chain and HKY simulation
  
  # SETUP
  dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
  set.seed(seed)
  
  options("logr.notes" = FALSE)
  log_open(file_name = paste0(out_dir, "simulation_log"))
  
  if(!is.null(seed)) {
    log_print(paste("Simulation seed set to:", seed)) 
  }
  
  transition_matrix <- as.matrix(
    read.csv("output/trade_matrix_daily_probabilities_eu.csv",
             row.names = 1,
             check.names = FALSE))
  
  # TRANSMISSION SIMULATION USING NOSOI
  trans_simulation <- run_nosoi(transition_matrix, max_infections = max_infections, sim_length = 365, out_dir = out_dir)
  
  # RUN HKY SUBSTITION
  host_data <-  getHostData(trans_simulation)
  ref_genome <- read.fasta("data/ref_genome.fasta", forceDNAtolower = FALSE, set.attributes = FALSE)[[1]]
  
  simulation_hky <- hky_nosoi(ref_genome = ref_genome, host_data = host_data)
  simulation_hky$seq <- sapply(simulation_hky$seq, function(x) paste(x, collapse = "")) # Collapse character vectors to strings
  
  # Save sequences
  simulation_hky$date.time <- simulation_hky$inf.time / 365.25 # Create a BEAST friendly time format
  IDs <- paste(simulation_hky$hosts.ID, simulation_hky$current.in, simulation_hky$date.time, sep = "|")
  sequences <- simulation_hky$seq
  names(sequences) <- IDs
  multifasta <- Biostrings::DNAStringSet(sequences)
  Biostrings::writeXStringSet(multifasta, paste0(out_dir, "sequences.fasta"))
  
  # Final output
  write.csv(simulation_hky, paste0(out_dir, "simulation_data.csv"))
  
  log_print("Simulation finished successfully. Final output saved as simulation_data.csv")
  log_close()
}

# Arguments (positional with defaults, no checks)
args <- commandArgs(trailingOnly = TRUE)
cat("Arguments received:", length(args), "->", paste(args, collapse = ", "), "\n")
max_infections <- if(length(args) >= 1) as.numeric(args[1])   else 10000
sim_length     <- if(length(args) >= 2) as.numeric(args[2])   else 365
seed           <- if(length(args) >= 3) as.numeric(args[3])   else NULL
out_dir        <- if(length(args) >= 4) as.character(args[4]) else "/output/simulation/"

#setwd("..")
base_dir <- normalizePath(".")
out_dir <- paste0(base_dir, out_dir)

cat("Working directory:", getwd(), "\n")
cat("Output directory:", out_dir, "\n")

seed = 12092 # For testing

run_simulation(max_infections = max_infections,
               sim_length = sim_length,
               seed = seed,
               out_dir = out_dir)

