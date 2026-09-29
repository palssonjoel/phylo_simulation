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
library(gridExtra)

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
  
  
  target <- setNames(pop_stats$herd_proportion, pop_stats$exporter)    
  
  stat <- function(M) { v <- Re(eigen(t(M))$vectors[,1]); setNames(v/sum(v), rownames(M)) }
  pi_stat <- stat(transition_matrix)
  
  ratio <- target / pi_stat                     # >1 = country is under-represented by trade
  lr <- log(ratio)
  
  w <- 1 + 0.15 * tanh(lr / 4)        
  R_base <- 0.95
  R_local <- R_base * w
  R_local <- round(R_local, 2)
  
  
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
  
  states <- colnames(transition_matrix)
  
  if (!identical(sort(states), sort(names(R_local)))) {
    stop(
      "State mismatch!\n",
      "Transition matrix states: ",
      paste(sort(states), collapse = ", "), "\n",
      "R_local states: ",
      paste(sort(names(R_local)), collapse = ", "), "\n",
      "Missing from R_local: ",
      paste(setdiff(states, names(R_local)), collapse = ", "), "\n",
      "Extra in R_local: ",
      paste(setdiff(names(R_local), states), collapse = ", ")
    )
  }
  
  if (anyNA(R_local)) {
    stop("R_local contains NA values.")
  }
  
  p_trans_func_diff <- function(t, current.in, p_max, t_incubation) {
    
    R <- NULL
    
    # Nosoi sanity checks demands that each state is explicitly written,
    # hence this nightmare >:(
    if (current.in == "Albania") R <- R_local[current.in]
    if (current.in == "Austria") R <- R_local[current.in]
    if (current.in == "Belgium") R <- R_local[current.in]
    if (current.in == "Bulgaria") R <- R_local[current.in]
    if (current.in == "Croatia") R <- R_local[current.in]
    if (current.in == "Cyprus") R <- R_local[current.in]
    if (current.in == "Czechia") R <- R_local[current.in]
    if (current.in == "Denmark") R <- R_local[current.in]
    if (current.in == "Estonia") R <- R_local[current.in]
    if (current.in == "Finland") R <- R_local[current.in]
    if (current.in == "France") R <- R_local[current.in]
    if (current.in == "Germany") R <- R_local[current.in]
    if (current.in == "Greece") R <- R_local[current.in]
    if (current.in == "Hungary") R <- R_local[current.in]
    if (current.in == "Ireland") R <- R_local[current.in]
    if (current.in == "Italy") R <- R_local[current.in]
    if (current.in == "Latvia") R <- R_local[current.in]
    if (current.in == "Lithuania") R <- R_local[current.in]
    if (current.in == "Luxembourg") R <- R_local[current.in]
    if (current.in == "Malta") R <- R_local[current.in]
    if (current.in == "Netherlands") R <- R_local[current.in]
    if (current.in == "Poland") R <- R_local[current.in]
    if (current.in == "Portugal") R <- R_local[current.in]
    if (current.in == "Romania") R <- R_local[current.in]
    if (current.in == "Serbia") R <- R_local[current.in]
    if (current.in == "Slovakia") R <- R_local[current.in]
    if (current.in == "Slovenia") R <- R_local[current.in]
    if (current.in == "Spain") R <- R_local[current.in]
    if (current.in == "Sweden") R <- R_local[current.in]
    if (current.in == "Switzerland") R <- R_local[current.in]
    
    if (is.null(R)) {
      stop("Unknown state: ", current.in)
    }
    
    if (t < t_incubation) {
      return(0)
    }
    
    return(p_max * R)
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
  t_incub <- t_incub_func(n)
  p_max <- p_max_func(n)
  n_contacts <- replicate(n, n_contact_func(0))
  p_moves <- replicate(n, p_move_func(0))
  
  # Histograms
  incub_plot <- ggplot(data.frame(t_incub), aes(x = t_incub)) +
    geom_histogram(bins = 50) +
    labs(
      title = "Distribution of incubation time",
      x = "Incubation time",
      y = "Frequency"
    ) +
    theme_minimal()
  
  p_max_plot <- ggplot(data.frame(p_max), aes(x = p_max)) +
    geom_histogram(bins = 50) +
    labs(
      title = "Distribution of maximum transmission probability",
      x = "p_max",
      y = "Frequency"
    ) +
    theme_minimal()
  
  n_contacts_plot <- ggplot(data.frame(n_contacts), aes(x = n_contacts)) +
    geom_histogram(
      breaks = seq(-0.5, max(n_contacts) + 0.5, by = 1)
    ) +
    labs(
      title = "Distribution of nContact",
      x = "Number of contacts",
      y = "Frequency"
    ) +
    theme_minimal()
  
  movement_plot <- ggplot(data.frame(p_moves), aes(x = p_moves)) +
    geom_histogram() +
    labs(
      title = "Distribution of pMove",
      x = "Move probability",
      y = "Frequency"
    ) +
    coord_cartesian(xlim = c(0, 0.1)) +
    theme_minimal()
  
  combined_plot <- grid.arrange(
    incub_plot,
    p_max_plot,
    n_contacts_plot,
    movement_plot,
    ncol = 2
  )
  
  ggsave(
    filename = paste0(out_dir,"histograms.png"),
    plot = combined_plot,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
  
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
                         
                         pTrans = p_trans_func_diff,
                         param.pTrans = list(p_max=p_max_func,t_incubation=t_incub_func),
                         timeDep.pTrans=FALSE,
                         diff.pTrans=TRUE,
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
    
    # Get and save transmission tree
    print("Getting transmission tree...")
    tree <- getTransmissionTree(trans_simulation)
    write.beast(tree, file.path(simulation_dir, "transmission_tree.nexus"))
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