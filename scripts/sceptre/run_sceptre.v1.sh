#!/bin/bash
#SBATCH --job-name=sceptre_v1
#SBATCH --partition=highmem_30d_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=700G
#SBATCH --time=29-23:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=hx37930@uga.edu
#SBATCH --mail-type=END,FAIL

# load R
ml R/4.3.1-foss-2022a

outputDir="/scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/merged_batch12/default/trans"
# file path of gene info
all_genes_fp="/scratch/hx37930/project/MTAG/06.crispr/results/scRNA_crispri_merged/combined_parents/yk3/DGE_filtered/all_genes.csv"
# file path of gene matrix
gene_mat_fp="/scratch/hx37930/project/MTAG/06.crispr/results/scRNA_crispri_merged/combined_parents/yk3/DGE_filtered/DGE.mtx"
# file path of gRNAs matrix
grna_mat_fp="/scratch/hx37930/project/MTAG/06.crispr/results/scRNA_crispri_merged/combined_crispr/yk3/guide_RNAs_filtered/count_matrix.mtx"
# directory of gRNA info
grna_dir="/scratch/hx37930/project/MTAG/06.crispr/results/scRNA_crispri_merged/combined_crispr/yk3/guide_RNAs_filtered"
all_grnas_fp_raw="${grna_dir}/all_guides.csv"
all_grnas_fp_clean="${grna_dir}/all_guides.clean.csv"
# file path of gRNA information, adding chromosome and pisition information
guide_chr="/scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/guideList.chr.uniq.txt"
# moi: high, or low
moi="high"
# remove cells that contain high proportioni (20%) of mitochrondria genes
p_mito_threshold=0.2
# run mode: cis, trans, both
mode="trans"  

if [ $mode = both ]; then
	mkdir -p ${outputDir}/cis ${outputDir}/trans
fi

# remove the 1st column and rename header
awk 'BEGIN{FS=OFS=","}NR==1{print "grna_id","grna_target","genome"}NR>1{print $2,$3,$4}' ${all_grnas_fp_raw} > ${all_grnas_fp_clean}

scriptDir="/scratch/hx37930/project/MTAG/scripts/scRNA_crispri"
Rscript ${scriptDir}/run_sceptre.v1.R --outputDir ${outputDir} --geneInfo ${all_genes_fp} --guideInfo ${all_grnas_fp_clean} --geneMatrix ${gene_mat_fp} --guideMatrix ${grna_mat_fp} --guideChr ${guide_chr} --moi ${moi} --mitoProp ${p_mito_threshold} --type ${mode}
