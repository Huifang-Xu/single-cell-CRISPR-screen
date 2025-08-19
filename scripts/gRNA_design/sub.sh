#!/bin/bash
#SBATCH --job-name=gRNA		# Job name (testBowtie2)
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=4	 	# CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=50G			# Memory per node (4GB); by default using M as unit
#SBATCH --time=1-23:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=NONE                   # Do not export any user’s explicit environment variables to compute node
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=hx37930@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

###############script##########
ml ANTLR/2.7.7-GCCcore-11.3.0-Java-11

# assign input and output
dataDir="/scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/finemAll_800bp_08112024"
reference="/scratch/hx37930/project/MTAG/06.crispr/reference/hg38.fa"
database="/scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/hg38_cas9ngg_database"
outDir="/scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/finemAll_800bp_08112024"
inFile="/scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/finemAll_800bp_08112024/fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.top3.tissue.400bp.redesignRNAs.txt"
output_prefix="fine_mappedSNP.PIP05.removeTwoRegions"

cd ${outDir}

# create a bed file to extract a 200 bp region centering each fine-mapped SNPs
# create a bed file
awk 'BEGIN{FS=OFS="\t"}NR>1{print $1":"strtonum($2)-400"-"strtonum($2)+400}' ${inFile} > ${dataDir}/${output_prefix}.bed
# extract target regions
ml SAMtools/1.17-GCC-12.2.0
samtools faidx $reference -r ${dataDir}/${output_prefix}.bed -o  ${dataDir}/${output_prefix}.fasta

# create off-target sites database
#java -Xmx4g -jar /home/hx37930/software/FlashFry/FlashFry-assembly-1.15.jar \
#        index \
#        --tmpLocation ${outDir}/tmp \
#        --database ${outDir}/hg38_cas9ngg_database \
#        --reference $reference \
#        --enzyme spcas9ngg

# discover candidate targets and their potential off-target effets
java -Xmx4g -jar /home/hx37930/software/FlashFry/FlashFry-assembly-1.15.jar \
	discover \
	--database ${database} \
	--fasta ${dataDir}/${output_prefix}.fasta \
	--output ${outDir}/${output_prefix}.gRNAs.txt

# score the discovered sites
java -Xmx4g -jar /home/hx37930/software/FlashFry/FlashFry-assembly-1.15.jar \
        score \
        --input ${outDir}/${output_prefix}.gRNAs.txt \
        --output ${outDir}/${output_prefix}.gRNAs.scored.txt \
        --scoringMetrics JostandSantos,hsu2013,doench2014ontarget,doench2016cfd,dangerous,minot \
        --database ${database}

