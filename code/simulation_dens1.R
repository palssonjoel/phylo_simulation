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
#library(gridExtra)

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
  
  # -------------------------------------------------------------------------
  # Country-specific transmission scaling
  #
  # The movement matrix defines the long-run distribution of animals across
  # countries (pi_stat). This is compared with the observed proportion of the
  # total EU herd in each country (target).
  #
  # target / pi_stat therefore measures whether a country is under- or
  # over-represented by the movement network relative to its herd size:
  #   >1  = larger herd share than implied by movement
  #   <1  = smaller herd share than implied by movement
  #
  # Transmission scaling is limited to approximately +/-15%.
  # -------------------------------------------------------------------------
  
  # Define target population distribution
  target <- setNames(
    pop_stats$herd_proportion,
    pop_stats$exporter
  )
  
  # Stationary distribution of the animal movement network
  stat <- function(M) {
    v <- Re(eigen(t(M))$vectors[, 1])
    setNames(v / sum(v), rownames(M))
  }
  
  pi_stat <- stat(transition_matrix)
  
  # Relative herd representation compared with movement-network representation
  ratio <- target / pi_stat
  
  # Log-ratio: 0 means herd share and movement share are equal
  lr <- log(ratio)
  
  # Bounded country-specific transmission multiplier
  w <- 1 + 0.15 * tanh(lr / 4)
  
  # Baseline transmission scaling and country-specific adjustment
  R_base <- 0.95
  R_local <- R_base * w
  
  # Round for reproducibility/readability
  R_local <- round(sort(R_local), 2)
  
  
  #############################################################################
  #                               DEFINE NOSOI FUNCTIONS
  #############################################################################
  
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
  
  # Beta distribution
  p_move_func <- function(t) {
    return(rbeta(1, shape1 = 1, shape2 = 99))
  }
  
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
  
  p_max_func <- function(x) {
    pmax(0, rbeta(x, shape1 = 3, shape2 = 10))
  }
  
  p_trans_func_diff <- function(t, current.in, p_max, t_incubation) {
    # This vairation of the p_trans function produces different p_trans values
    # for different demes, based on the difference between the distribution
    # of animals inferred from the transition matrix, and the observed 
    # distribution in the data.
    
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
  
  p_trans_density1 <- function(t, current.in, host.count, p_max, t_incubation) {
    
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
    
    if (t < t_incubation) {
      return(0)
    }
    
    # Introduce a density dependent suppression factor
    K <- 1000   # Scale at which suppression starts
    h <- 1      # Sharpness of suppression
    
    density_factor <- 1 / (1 + (host.count / K)^h)
    
    p_trans <- p_max * R * density_factor
    
    return(p_trans)
  }
  
  p_trans_density2 <- function(t, current.in, host.count, p_max, t_incubation) {
    # This adds a density dependent suppression to p_trans_func_diff(),
    # which reduces the number of active hosts in each country after 
    # increasing above the cutoff point K. This is introduced to 
    # allow for endemic spread, without the epidemic scaling to unsustainable
    # levels when expanded to run over 10-20 years.
    
    if (current.in == "Albania") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Austria") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Belgium") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Bulgaria") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Croatia") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Cyprus") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Czechia") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Denmark") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Estonia") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Finland") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "France") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Germany") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Greece") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Hungary") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Ireland") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Italy") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Latvia") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Lithuania") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Luxembourg") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Malta") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Netherlands")pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Poland") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Portugal") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Romania") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Serbia") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Slovakia") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Slovenia") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Spain") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Sweden") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    if (current.in == "Switzerland") pop <- pop_stats$exporter_herd_size[pop_stats$exporter==current.in]
    
    if (t < t_incubation) {
      return(0)
    }
    
    # Density dependent suppression
    h <- 2   # sharpess of suppression. 1 = gradual, 2 = strong
    q <- 0.01 # prop of population to use for suppression
    K <- q * pop
    
    # Testing K = 1%: https://pmc.ncbi.nlm.nih.gov/articles/PMC4146608/?utm_source=chatgpt.com
    
    density_factor <- 1 / (1 + (host.count / K)^h)
    
    p_trans <- p_max * density_factor
    
    return(p_trans)
  }
  
  p_exit_func <- function(t, t_incubation) {
    if (t < t_incubation) { return(0) }
    else {
      return(1/7)  # ≈ 0.20/day → mean ~5 days post-incubation illness
    }
  }
  
  ################################################################################
  # Epidemic dynamics histograms
  ################################################################################
  # Mainly used when setting up, but not useful for actual runs. 
  
  # Number of simulated values
  n <- 10000
  
  # Generate values
  #t_incub <- t_incub_func(n)
  p_max <- p_max_func(n)
  #n_contacts <- replicate(n, n_contact_func(0))
  #p_moves <- replicate(n, p_move_func(0))
  
  # Histograms
  #incub_plot <- ggplot(data.frame(t_incub), aes(x = t_incub)) +
  #  geom_histogram(bins = 50) +
  #  labs(
  #    title = "Distribution of incubation time",
  #    x = "Incubation time",
  #    y = "Frequency"
  #  ) +
  #  theme_minimal()
  
  p_max_plot <- ggplot(data.frame(p_max), aes(x = p_max)) +
    geom_histogram(bins = 50) +
    labs(
      title = "Distribution of maximum transmission probability",
      x = "p_max",
      y = "Frequency"
    ) +
    theme_minimal()
  
  #n_contacts_plot <- ggplot(data.frame(n_contacts), aes(x = n_contacts)) +
  #  geom_histogram(
  #    breaks = seq(-0.5, max(n_contacts) + 0.5, by = 1)
  #  ) +
  #  labs(
  #    title = "Distribution of nContact",
  #    x = "Number of contacts",
  #    y = "Frequency"
  #  ) +
  #  theme_minimal()
  
  #movement_plot <- ggplot(data.frame(p_moves), aes(x = p_moves)) +
  #  geom_histogram() +
  #  labs(
  #    title = "Distribution of pMove",
  #    x = "Move probability",
  #    y = "Frequency"
  #  ) +
  #  coord_cartesian(xlim = c(0, 0.1)) +
   # theme_minimal()
  
  #combined_plot <- grid.arrange(
  #  incub_plot,
  #  p_max_plot,
  #  n_contacts_plot,
  #  movement_plot,
  #  ncol = 2
  #)
  
  #ggsave(
  #  filename = paste0(out_dir,"histograms.png"),
  #  plot = combined_plot,
  #  width = 10,
  #  height = 8,
  #  units = "in",
  #  dpi = 300
  #)
  
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
                         
                         pTrans = p_trans_density1,
                         param.pTrans = list(p_max=p_max_func,t_incubation=t_incub_func),
                         timeDep.pTrans=FALSE,
                         diff.pTrans=TRUE,
                         hostCount.pTrans = TRUE,
                         
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
# HKY model (NOT USED IN THIS VERSION)
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
  
  # SETUP
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  options("logr.notes" = FALSE)
  log_open(file_name = paste0(out_dir, "simulation"))
  
  transition_matrix <- as.matrix(
    read.csv("../output/trade_matrix_daily_probabilities_eu.csv",
             row.names = 1,
             check.names = FALSE))
  
  successful_runs <- 0
  attempts <- 0
  
  max_attempts <- run_no * 100
  
  # Run simulation until run_no of complete simulations are done, 
  # discarding failed runs. 
  while (successful_runs < run_no && attempts < max_attempts) {
    
    attempts <- attempts + 1
    
    seed <- sample.int(999999999, 1)
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
      
      #message(msg)
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
    final_time <- max(host_table$out.time, na.rm = TRUE)
    
    # Warn if the simulation reaches the host limit, but do not reject it
    if (n_hosts >= max_infections) {
      msg <- paste0(
        "Attempt ", attempts,
        " reached max host threshold (", max_infections,
        ") at day ", final_time, "."
      )
      log_print(msg)
    }
    
    # Accept simulations that reach at least 80% of the requested time
    min_acceptable_time <- 0.80 * sim_length
    if (final_time < min_acceptable_time && n_hosts <= 0.9 * max_infections) {
      msg <- paste0(
        "Attempt ", attempts,
        " rejected: epidemic ended at day ",
        round(final_time, 1),
        " (< 90% of ", sim_length, " days = ",
        round(min_acceptable_time, 1), " days)."
      )
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
    
    #message(msg)
    log_print(msg)
    
    # --------------------------------------------------
    # Save successful nosoi output
    # --------------------------------------------------
    
    saveRDS(
      trans_simulation,
      paste0(simulation_dir, "/nosoi_sim.rds")
    )
    print("Simulation saved!")
    
    # Get and save transmission tree
    start_time <- Sys.time()
    print("Getting transmission tree...")
    tree <- getTransmissionTree(trans_simulation)
    write.beast(tree, file.path(simulation_dir, "transmission_tree.nexus"))
    #write.beast.newick(tree, file.path(out_dir, "transmission_tree.nwk"))
    end_time <-  Sys.time()
    print("Transmission tree done! Saved as .nexus file")
    duration <- difftime(end_time, start_time, units = "mins")
    print(paste0("Time take: ", duration))
    
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

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

cat("Working directory:", getwd(), "\n")
cat("Output directory:", out_dir, "\n")

run_simulation(max_infections = max_infections,
               sim_length = sim_length,
               run_no = run_no,
               out_dir = out_dir)

