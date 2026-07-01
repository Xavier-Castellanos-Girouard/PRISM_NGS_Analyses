#!/bin/bash
#SBATCH --job-name=NestedPCR_v3          # Job name
#SBATCH --output=NestedPCR_v3_%j.out     # Standard output log (%j = job ID)
#SBATCH --error=NestedPCR_v3_%j.err      # Standard
#SBATCH --time=4:00:00                   # Max runtime (hh:mm:ss)
#SBATCH --cpus-per-task=12                # Number of CPU cores
#SBATCH --mem-per-cpu=16G

# Load required modules (if your HPC uses module system)
#module load bowtie2/2.3.5.1
#module load samtools/1.9
#module load python/3.10
#module load cutadapt/5.0

#
source ~/miniconda3/etc/profile.d/conda.sh

# Activate your conda environment (if needed)
conda activate bio_env

# Run your nested PCR pipeline
bash NestedPCR_MasterScript.sh
