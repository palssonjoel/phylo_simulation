#!/bin/bash
#SBATCH --job-name=seq_evolution_08102026
#SBATCH --partition=project
#SBATCH --output=seq_evolution_%j.out
#SBATCH --error=seq_evolution_%j.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

source /dpssi/data/Tools/Conda/latest/etc/profile.d/conda.sh
conda activate phylo

iqtree3 --alisim simulated_alignment \
    --root-seq ../config/Aswine_Arnsberg_6554_1979_merged.fasta,'A/swine/Arnsberg/6554/1979|EPI_ISL_539869' \
    -p ../config/partitions.nex \
    -t ../output/reduced_tree.nwk \
    --branch-scale 2.0808e-5
