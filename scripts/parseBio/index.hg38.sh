#!/bin/bash
#SBATCH --job-name=index          # Job name (testBowtie2)
#SBATCH --partition=batch               # Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1                      # Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=1               # CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=50G                        # Memory per node (4GB); by default using M as unit
#SBATCH --time=2-00:00:00                  # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=NONE                   # Do not export any user’s explicit environment variables to compute node
#SBATCH --output=%x_%j.out              # Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err               # Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=hx37930@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL            # Mail events (BEGIN, END, FAIL, ALL)

# activate conda env
source /home/hx37930/Miniconda3/etc/profile.d/conda.sh
conda activate spipe

workDir="/scratch/hx37930/reference/hg38"

split-pipe \
	--mode mkref \
	--genome_name hg38 \
	--fasta ${workDir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
	--genes ${workDir}/Homo_sapiens.GRCh38.109.gtf.gz \
	--output_dir ${workDir} --start_timeout 0
