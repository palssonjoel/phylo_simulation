#!/bin/bash

#SBATCH --partition=project
sbatch phylo_simulation_1y.slurm
sbatch phylo_simulation_5y.slurm
sbatch phylo_simulation_10y.slurm
sbatch phylo_simulation_20y.slurm
