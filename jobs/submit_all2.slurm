#!/bin/bash

#SBATCH --partition=project
sbatch phylo_simulation_1y_2.slurm
sbatch phylo_simulation_5y_2.slurm
sbatch phylo_simulation_10y_2.slurm
sbatch phylo_simulation_20y_2.slurm
