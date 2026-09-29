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
  
  # Deme pig population stats for density dependence 
  pop_stats <- read.csv("../output/trade_movement_rates_eu.csv") |> 
    select(exporter, exporter_herd_size) |> 
    distinct()
  
  # Total population and deme proportion of total
  pig_total <- sum(pop_stats$exporter_herd_size)
  pop_stats$herd_proportion <- pop_stats$exporter_herd_size / pig_total
  
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
  
  # pMove is probability that one pig moves between locations. 
  # Beta distribution
  p_move_func <- function(t) {
    return(rbeta(1, shape1 = 1, shape2 = 99))
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
    abs(round(rnorm(1, 1, 1)))
  }
  
  # pTrans depends on incubation time, which is often cited between 1-3 days
  # p_max is a constant probability of transmission
  p_trans_func <- function(t, p_max, t_incubation) {
    if(t < t_incubation){p = 0}
    if(t >= t_incubation){p = p_max}
    return(p)
  }
  
  t_incub_func <- function(x){pmax(0, rnorm(x, mean = 2, sd = 0.5))}
  #p_max_func <- function(x){rbeta(x, shape1=1, shape2=3)}
  
  # Mean = 0.167
  p_max_func <- function(x) {
    pmax(0.05, rbeta(x, shape1 = 1, shape2 = 10))
  }
  
  p_trans_func_ <- function(t, current.in, host.count) {
    
    if (current.in == "Albania") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Albania"]
    } else if (current.in == "Austria") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Austria"]
    } else if (current.in == "Belgium") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Belgium"]
    } else if (current.in == "Bulgaria") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Bulgaria"]
    } else if (current.in == "Croatia") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Croatia"]
    } else if (current.in == "Cyprus") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Cyprus"]
    } else if (current.in == "Czechia") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Czechia"]
    } else if (current.in == "Denmark") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Denmark"]
    } else if (current.in == "Estonia") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Estonia"]
    } else if (current.in == "Finland") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Finland"]
    } else if (current.in == "France") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "France"]
    } else if (current.in == "Germany") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Germany"]
    } else if (current.in == "Greece") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Greece"]
    } else if (current.in == "Hungary") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Hungary"]
    } else if (current.in == "Ireland") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Ireland"]
    } else if (current.in == "Italy") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Italy"]
    } else if (current.in == "Latvia") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Latvia"]
    } else if (current.in == "Lithuania") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Lithuania"]
    } else if (current.in == "Luxembourg") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Luxembourg"]
    } else if (current.in == "Malta") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Malta"]
    } else if (current.in == "Netherlands") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Netherlands"]
    } else if (current.in == "Poland") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Poland"]
    } else if (current.in == "Portugal") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Portugal"]
    } else if (current.in == "Romania") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Romania"]
    } else if (current.in == "Serbia") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Serbia"]
    } else if (current.in == "Slovakia") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Slovakia"]
    } else if (current.in == "Slovenia") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Slovenia"]
    } else if (current.in == "Spain") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Spain"]
    } else if (current.in == "Sweden") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Sweden"]
    } else if (current.in == "Switzerland") {
      herd_prop <- pop_stats$herd_proportion[pop_stats$exporter == "Switzerland"]
    } else {
      stop("Unknown state: ", current.in)
    }
    # Deme-specific carrying capacity
    K_deme <- max_infections * herd_prop
    # Density-dependent suppression
    sup_factor <- max(0.1, 1 - host.count / K_deme)
    # Incubation period
    if (t < 2) {
      return(0)
    }
    # Transmission probability
    return(0.5 * sup_factor)
  }
  
  
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
      return(1/10)  # ≈ 0.20/day → mean ~5 days post-incubation illness
    }
  }
  
  ################################################################################
  # Epidemic dynamics histograms
  ################################################################################
  # Mainly used when setting up, but not useful for actual runs. 
  # Remove comments to save them for each iteration.
  
  # Number of simulated values
  n <- 10000
  
  # Generate values
  #t_incub <- t_incub_func(n)
  p_max <- p_max_func(n)
  #n_contacts <- replicate(n, n_contact_func(0))
  #p_moves <- replicate(n, p_move_func(0))
  
  # Histograms
 #incub_hist <- hist(
#    t_incub,
#    breaks = 50,
#    main = "Distribution of incubation time",
#    xlab = "Incubation time",
#    ylab = "Frequency"
#  )
  
 p_max_hist <- hist(
    p_max,
    breaks = 50,
    main = "Distribution of maximum transmission probability",
    xlab = "p_max",
    ylab = "Frequency"
  )
 
 
  
  #n_contacts_hist <- hist(
  #  n_contacts,
  #  breaks = seq(-0.5, max(n_contacts) + 0.5, by = 1),
  #  main = "Distribution of nContact",
  #  xlab = "Number of contacts",
  #  ylab = "Frequency"
  #)
  
  #movement_hist <- hist(
  #  p_moves,
  #  #breaks = 50,
  #  main = "Distribution of pMove",
  #  xlab = "Move probability",
  #  ylab = "Frequency",
  #  xlim = c(0,0.1)
  #)
  
  
  # Save plots
  #png(file.path(out_dir, "incub_hist.png"), width = 800, height = 600)
  #plot(incub_hist, main = "Distribution of incubation time", xlab = "Incubation time")
  #dev.off()
  
  #png(file.path(out_dir, "p_max_hist.png"), width = 800, height = 600)
  #plot(p_max_hist, main = "Distribution of maximum transmission probability", xlab = "p_max")
  #dev.off()
  
  #png(file.path(out_dir, "n_contacts_hist.png"), width = 800, height = 600)
  #plot(n_contacts_hist, main = "Distribution of nContact", xlab = "Number of contacts")
  #dev.off()
  
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
                         hostCount.pTrans = FALSE,
                         
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

create_alignment <- function(hky_output, filename, out_dir) {
  # Extract sequences from the hky + nosoi output, and write to an aligned
  # fasta file. 
  
  simulation_hky <- hky_output
  simulation_dir <- out_dir
  
  simulation_hky$date.time <- simulation_hky$out.time / 365.25 # Create a BEAST friendly time format
  IDs <- paste(simulation_hky$hosts.ID, simulation_hky$current.in, simulation_hky$date.time, sep = "|")
  sequences <- simulation_hky$seq
  names(sequences) <- IDs
  multifasta <- Biostrings::DNAStringSet(sequences)
  Biostrings::writeXStringSet(multifasta, paste0(simulation_dir, "/", filename, ".fasta"))
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
                           run_no, 
                           out_dir) {
  # This runs both the transmission chain and HKY simulation
  
  # SETUP
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  #seed = 1242 
  #set.seed(seed) # For testing
  
  options("logr.notes" = FALSE)
  log_open(file_name = paste0(out_dir, "simulation"))
  
  transition_matrix <- as.matrix(
    read.csv("../output/trade_matrix_daily_probabilities_eu.csv",
             row.names = 1,
             check.names = FALSE))
  
  successful_runs <- 0
  attempts <- 0
  
  max_attempts <- run_no * 1000
  
  while (successful_runs < run_no && attempts < max_attempts) {
    
    attempts <- attempts + 1
    
    seed <- sample.int(9999999, 1)
    set.seed(seed)
    
    log_print(
      paste0(
        "Attempt ", attempts,
        " | Successful simulations: ",
        successful_runs, "/", run_no,
        " | Seed = ", seed
      )
    )
    
    simulation_dir <- paste0(
      out_dir, "/",
      successful_runs + 1
    )
    
    dir.create(
      simulation_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    # --------------------------------------------------
    # Run epidemic simulation
    # --------------------------------------------------
    
    trans_simulation <- tryCatch({
      
      run_nosoi(
        transition_matrix,
        max_infections = max_infections,
        sim_length = sim_length,
        out_dir = simulation_dir
      )
      
    }, error = function(e) {
      
      msg <- paste0(
        "Attempt ", attempts,
        " failed with error: ",
        conditionMessage(e)
      )
      
      message(msg)
      log_print(msg)
      
      return(NULL)
    })
    
    if (is.null(trans_simulation)) {
      next
    }
    
    # --------------------------------------------------
    # Check whether epidemic was successful
    # --------------------------------------------------
    
    host_table <- getTableHosts(trans_simulation)
    
    n_hosts <- nrow(host_table)
    
    final_time <- max(
      host_table$out.time,
      na.rm = TRUE
    )
    
    if (final_time < sim_length) {
      
      msg <- paste0(
        "Attempt ", attempts,
        " rejected: epidemic ended at day ",
        final_time,
        " (< ", sim_length, " days)."
      )
      
      message(msg)
      log_print(msg)
      
      next
    }
    
    # --------------------------------------------------
    # Successful simulation
    # --------------------------------------------------
    
    successful_runs <- successful_runs + 1
    
    msg <- paste0(
      "Successful simulation ",
      successful_runs,
      "/", run_no,
      " obtained after ",
      attempts,
      " attempts."
    )
    
    message(msg)
    log_print(msg)
    
    # --------------------------------------------------
    # Save successful nosoi output
    # --------------------------------------------------
    
    saveRDS(
      trans_simulation,
      paste0(simulation_dir, "/nosoi_sim.rds")
    )
    
    # --------------------------------------------------
    # Run HKY
    # --------------------------------------------------
    
    host_data <- getHostData(trans_simulation)
    
    ref_genome <- read.fasta(
      "../data/ref_genome.fasta",
      forceDNAtolower = FALSE,
      set.attributes = FALSE
    )[[1]]
    
    simulation_hky <- hky_nosoi(
      ref_genome = ref_genome,
      host_data = host_data
    )
    
    simulation_hky$seq <- sapply(
      simulation_hky$seq,
      function(x) paste(x, collapse = "")
    )
    
    create_alignment(
      simulation_hky,
      "full_seqs",
      simulation_dir
    )
    
    write.csv(
      simulation_hky,
      paste0(
        simulation_dir,
        "/simulation_data.csv"
      )
    )
    
    # Get and save transmission tree
    print("Getting transmission tree...")
    tree <- getTransmissionTree(trans_simulation)
    write.beast(tree, file.path(out_dir, "transmission_tree.nexus"))
    #write.beast.newick(tree, file.path(out_dir, "transmission_tree.nwk"))
    print("Transmission tree done! Saved as .nexus file")
    
    log_print(
      paste0(
        "Simulation ",
        successful_runs,
        " complete."
      )
    )
    
    log_print(
      "========================================================================="
    )
  }
  
  # --------------------------------------------------
  # Final status
  # --------------------------------------------------
  
  if (successful_runs < run_no) {
    
    warning(
      "Only ",
      successful_runs,
      " successful simulations obtained after ",
      attempts,
      " attempts."
    )
    
  } else {
    
    log_print(
      paste0(
        "All ",
        run_no,
        " successful simulations completed after ",
        attempts,
        " attempts."
      )
    )
  }
}
  
# Arguments (positional with defaults, no checks)
args <- commandArgs(trailingOnly = TRUE)
cat("Arguments received:", length(args), "->", paste(args, collapse = ", "), "\n")
max_infections <- if(length(args) >= 1) as.numeric(args[1])   else 10000
sim_length     <- if(length(args) >= 2) as.numeric(args[2])   else 365
run_no         <- if(length(args) >= 3) as.numeric(args[3])   else 50
out_dir        <- if(length(args) >= 4) as.character(args[4]) else "/output/simulation/"

#setwd("..")
base_dir <- normalizePath("..")
out_dir <- paste0(base_dir, out_dir)

cat("Working directory:", getwd(), "\n")
cat("Output directory:", out_dir, "\n")

run_simulation(max_infections = max_infections,
               sim_length = sim_length,
               run_no = run_no,
               out_dir = out_dir)

