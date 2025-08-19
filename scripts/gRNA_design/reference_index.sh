#!/bin/bash
#SBATCH --job-name=index
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=6-24:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=hx37930@uga.edu
#SBATCH --mail-type=END,FAIL

#ml Bowtie/1.3.1-GCC-10.2.0

#bowtie /scratch/hx37930/reference/hg19.fa /scratch/hx37930/reference


ml SAMtools/1.16.1-GCC-11.3.0
samtools faidx /scratch/hx37930/project/MTAG/06.crispr/reference/hg38.fa
samtools faidx /scratch/hx37930/project/MTAG/06.crispr/reference/hg38.fa.masked
