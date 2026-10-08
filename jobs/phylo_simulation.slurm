#!/bin/bash
#SBATCH --job-name=phylo_simulation_10000_1096
#SBATCH --partition=project
#SBATCH --output=/dpssi/data/Projects/inf_fluzoomark/proj/JOELP/phylo_simulation/output/job_logs/phylo_simulation_10000_1096.out
#SBATCH --error=/dpssi/data/Projects/inf_fluzoomark/proj/JOELP/phylo_simulation/output/job_logs/phylo_simulation_10000_1096.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

# Load conda
source /dpssi/data/Tools/Conda/latest/etc/profile.d/conda.sh
conda activate R

# Go to project directory
cd /dpssi/data/Projects/inf_fluzoomark/proj/JOELP/phylo_simulation/code || exit 1

# Run simulation
Rscript simulation_only.R 10000 1096 1 ../output/simulation_10000_1096/
