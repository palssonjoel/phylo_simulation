#!/bin/bash
#SBATCH --job-name=seq_evolution
#SBATCH --partition=project
#SBATCH --output=seq_evolution_07102026.out
#SBATCH --error=seq_evolution_07102026.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

source /dpssi/data/Tools/Conda/latest/etc/profile.d/conda.sh
conda activate phylo

iqtree3 --alisim Aswine_Arnsberg_6554_1979_merged.fasta \
    -p partitions.nex \
    -t reduced_tree.nwk \
    --branch-scale 2.0808e-5
