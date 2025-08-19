#!/bin/bash
#SBATCH --job-name=comb_libs          # Job name (testBowtie2)
#SBATCH --partition=batch              # Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1                      # 1 task (process) for below commands
#SBATCH --cpus-per-task=10               # CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=50G                        # Memory per node (4GB); by default using M as unit
#SBATCH --time=1-23:00:00                  # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out              # Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=hx37930@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL                 # Mail events (BEGIN, END, FAIL, ALL)

#load the conda environment
source /home/hx37930/Miniconda3/etc/profile.d/conda.sh
conda activate spipe

sublibs_parent_list="/scratch/hx37930/project/MTAG/scripts/scRNA_crispri/sublibs.parent.list"
sublibs_crispr_list="/scratch/hx37930/project/MTAG/scripts/scRNA_crispri/sublibs.crispr.list"
combine_parent_dir="/scratch/hx37930/project/MTAG/06.crispr/results/scRNA_crispri_merged_dCas9/combined_parents"
combine_crispr_dir="/scratch/hx37930/project/MTAG/06.crispr/results/scRNA_crispri_merged_dCas9/combined_crispr"


#---------
#Set which
#steps run
#---------
step3=true
#---------

#######################################################
########## STEP 3. combine all sublibraries ###########
#######################################################
if [ $step3 = true ]; then

echo "-=-=-=-=-=-=-=-STEP 3-=-=-=-=-=-=-=-\n\n"
# combine all sublibraries without parent
split-pipe \
    --mode comb \
    --sublib_list ${sublibs_parent_list} \
    --output_dir ${combine_parent_dir}

# combine all sublibraries with a parent
split-pipe \
    --mode comb \
    --parent_dir ${combine_parent_dir} \
    --output_dir ${combine_crispr_dir} \
    --sublib_list ${sublibs_crispr_list}
fi

