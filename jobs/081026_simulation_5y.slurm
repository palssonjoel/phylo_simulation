#!/bin/bash
#SBATCH --job-name=081026_sim_5y
#SBATCH --partition=project
#SBATCH --output=081026_sim_5y_%j.out
#SBATCH --error=081026_sim_5y_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

source /dpssi/data/Tools/Conda/latest/etc/profile.d/conda.sh
conda activate R

cd /dpssi/data/Projects/inf_fluzoomark/proj/JOELP/phylo_simulation/code || exit 1

mkdir -p ../output/simulation_5y

Rscript simulation_dens2.R 100000000 1825 1 ../output/simulation_5y/
