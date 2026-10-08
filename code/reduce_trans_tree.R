library(nosoi)
library(dplyr)
library(ggplot2)
library(patchwork)
library(treeio)


reduce_tree <- function(sim, out_dir) {
  
  seed <- sample.int(999999999, 1)
  set.seed(seed)
  print(paste0("Using seed: ", seed))
  
  cat("Working directory:", getwd(), "\n")
  cat("Output directory:", out_dir, "\n")
  
  host_table <- getTableHosts(sim)
  host_table$date.time <- host_table$out.time / 365.25 # Crete beast-friendly time format
  
  # Reduce host_table, include first and last cases
  index <- slice_head(host_table)
  final <- slice_tail(host_table)
  
  # Random sampling
  host_subsample <- slice_sample(host_table, n=10000); host_subsample <- rbind(index, host_subsample, final)
  
  # Break time into bins and sample proportionally in each bin
  #n_sample <- 100
  
  #host_subsample <- host_table %>%
  #  mutate(time_bin = cut(inf.time, breaks = 10)) %>%
  #  group_by(time_bin, inf.in) %>%
  #  slice_sample(prop = 0.01) %>%
  #  ungroup()
  
  # 1. sample hosts then add all ancestors
  ht <- as.data.table(host_table); setkey(ht, hosts.ID)
  
  keep <- unique(host_subsample$hosts.ID); frontier <- keep
  repeat {
    par <- unique(ht[.(frontier), inf.by])
    par <- par[!is.na(par) & par != "NA"]
    new <- setdiff(par, keep)
    if (!length(new)) break
    keep <- c(keep, new); frontier <- new
  }
  
  # 2-3. filter host table and state table to the same IDs
  sim_sub <- sim
  sim_sub$host.info.A$table.hosts <- sim$host.info.A$table.hosts[hosts.ID %chin% keep]
  sim_sub$host.info.A$table.state <- sim$host.info.A$table.state[hosts.ID %chin% keep]
  sim_sub$host.info.A$N.infected  <- length(keep)
  
  # Counts
  reduced_n <- sim_sub$host.info.A$table.hosts |> nrow()
  full_n    <- sim$host.info.A$table.hosts |> nrow()
  
  print(paste0("Full n hosts: ", full_n))
  print(paste0("Reduced n hosts: ", reduced_n))
  # Check dynamics for consistency
  #reduced_dyn <- getDynamic(sim_sub)
  #full_dyn    <- getDynamic(sim)
  
  # Active hosts per timestep, by state
  #reduced_dyn |> 
  #  group_by(state, t) |> 
  #  ggplot(aes(t, Count, color=state)) +
  #  geom_path(linewidth = 1) +
  #  labs(
  #    title = "Epidemic dynamics",
  #   y = "No. Active hosts",
  #    x = "time"
  #  ) 
  
  #full_dyn |> 
  #  group_by(state, t) |> 
  #  ggplot(aes(t, Count, color=state)) +
  #  geom_path(linewidth = 1) +
  #  labs(
  #    title = "Epidemic dynamics",
  #    y = "No. Active hosts",
  #    x = "time"
  #  ) 
  
  
  # 4. nosoi's own functions, now on a small object
  print("Getting transmission tree...")
  reduced_tree <- getTransmissionTree(sim_sub)
  print("Done!")
  
  # 5. Save files
  print("Saving files...")
  
  saveRDS(sim_sub, paste0(out_dir, "reduced_nosoi_sim.rds"))
  write.csv(sim_sub$host.info.A$table.hosts, paste0(out_dir, "reduced_host_table.csv"))
  write.beast.newick(reduced_tree, paste0(out_dir, "reduced_tree.nwk"))
  
  print(paste0("All files saved in directory: ", out_dir))
}

# Arguments (positional with defaults, no checks)
args <- commandArgs(trailingOnly = TRUE)
cat("Arguments received:", length(args), "->", paste(args, collapse = ", "), "\n")
simulation     <- if(length(args) >= 1) as.numeric(args[1])  
out_dir        <- if(length(args) >= 2) as.character(args[2]) 

if(is.null(out_dir)) {
  stop("No output directory given.")
}

if(is.null(simulation)) {
  stop("No simulation provided.")
}

#setwd("..")
out_dir <- normalizePath(out_dir, mustWork = TRUE)

cat("Working directory:", getwd(), "\n")
cat("Output directory:", out_dir, "\n")

reduce_tree(sim=simulation, out_dir=out_dir)