#!/bin/bash
#SBATCH --job-name=sub6          # Job name (testBowtie2)
#SBATCH --partition=batch              # Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1                      # 1 task (process) for below commands
#SBATCH --cpus-per-task=10               # CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=150G                        # Memory per node (4GB); by default using M as unit
#SBATCH --time=3-23:00:00                  # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out              # Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=hx37930@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL                 # Mail events (BEGIN, END, FAIL, ALL)

#load the conda environment
source /home/hx37930/Miniconda3/etc/profile.d/conda.sh
conda activate spipe

ml FastQC/0.11.9-Java-11

fq1_parent_batch1="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/batch1/19117FL-35-01-06_S85_L004_R1_001.fastq.gz"
fq2_parent_batch1="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/batch1/19117FL-35-01-06_S85_L004_R2_001.fastq.gz"
fq1_parent_batch2="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/batch2/19117FL-35-01-06_S6703_L001_R1_001.fastq.gz"
fq2_parent_batch2="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/batch2/19117FL-35-01-06_S6703_L001_R2_001.fastq.gz"
fq1_guide_batch1="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/batch1/19117FL-35-02-06_S93_L004_R1_001.fastq.gz"
fq2_guide_batch1="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/batch1/19117FL-35-02-06_S93_L004_R2_001.fastq.gz"
fq1_guide_batch2="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/batch2/19117FL-35-02-06_S6703_L001_R1_001.fastq.gz"
fq2_guide_batch2="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/batch2/19117FL-35-02-06_S6703_L001_R2_001.fastq.gz"
parent_fq1="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/merge/19117FL-35-01-06_R1_001.fastq.gz"
parent_fq2="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/merge/19117FL-35-01-06_R2_001.fastq.gz"
crispr_fq1="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/merge/19117FL-35-02-06_R1_001.fastq.gz"
crispr_fq2="/scratch/hx37930/project/MTAG/01.data/scCRISPRscreens/merge/19117FL-35-02-06_R2_001.fastq.gz"
parent_output_dir="/scratch/hx37930/project/MTAG/06.crispr/results/scRNA_crispri_merged_dCas9/sub6"
crispr_output_dir="/scratch/hx37930/project/MTAG/06.crispr/results/scRNA_crispri_merged_dCas9/sub6_crispr"
parfile="/scratch/hx37930/project/MTAG/scripts/scRNA_crispri/parameterFile.txt"
genomeDir="/scratch/hx37930/reference/hg38_dCas9"
guidesFile="/scratch/hx37930/project/MTAG/scripts/scRNA_crispri/guideList.uniq.csv"

#---------
#Set which
#steps run
#---------
step1=false
step2=true
step3=true
#---------

#######################################################
### STEP 1. merge fastq files for each sublibrary ###
#######################################################
if [ $step1 = true ]; then
echo "-=-=-=-=-=-=-=-STEP 1-=-=-=-=-=-=-=-\n\n"
cat ${fq1_parent_batch1} ${fq1_parent_batch2} > ${parent_fq1}
cat ${fq2_parent_batch1} ${fq2_parent_batch2} > ${parent_fq2}
cat ${fq1_guide_batch1} ${fq1_guide_batch2} > ${crispr_fq1}
cat ${fq2_guide_batch1} ${fq2_guide_batch2} > ${crispr_fq2}
fi

#######################################################
## STEP 2. process transcriptome for one sublibrary ###
#######################################################
if [ $step2 = true ]; then
echo "-=-=-=-=-=-=-=-STEP 2-=-=-=-=-=-=-=-\n\n"
split-pipe \
    --mode all \
    --parfile ${parfile} \
    --chemistry v2 \
    --genome_dir ${genomeDir} \
    --fq1 ${parent_fq1} \
    --output_dir ${parent_output_dir} \
    --sample yk3 A1-D12
fi

#######################################################
###### STEP 3. process crispr for one sublibrary ######
#######################################################
if [ $step3 = true ]; then

echo "-=-=-=-=-=-=-=-STEP 3-=-=-=-=-=-=-=-\n\n"
split-pipe \
    --mode all \
    --parfile ${parfile} \
    --chemistry v2 \
    --crispr \
    --crsp_guides ${guidesFile} \
    --parent_dir ${parent_output_dir} \
    --output_dir ${crispr_output_dir} \
    --fq1 ${crispr_fq1}
fi

