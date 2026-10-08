#!/bin/bash
#SBATCH --job-name=phylo_5y_dens2
#SBATCH --partition=project
#SBATCH --output=/dpssi/data/Projects/inf_fluzoomark/proj/JOELP/phylo_simulation/output/job_logs/phylo_5y_dens_%j.out
#SBATCH --error=/dpssi/data/Projects/inf_fluzoomark/proj/JOELP/phylo_simulation/output/job_logs/phylo_5y_dens_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

source /dpssi/data/Tools/Conda/latest/etc/profile.d/conda.sh
conda activate R

cd /dpssi/data/Projects/inf_fluzoomark/proj/JOELP/phylo_simulation/code || exit 1

mkdir -p ../output/dens2/simulation_5y_dens

Rscript simulation_dens2.R 10000000 1825 1 ../output/dens2/simulation_5y_dens/
