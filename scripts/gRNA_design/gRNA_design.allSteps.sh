#!bin/bash

# this script is used to design gRNAs using FlashFry

############################################################
# step 0: install FlashFry and run test data
############################################################
# ml java (tutorial used version 8)
ml ANTLR/2.7.7-GCCcore-11.3.0-Java-11
# download FlashFry
wget https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar

# download sample data: chr22
wget https://raw.githubusercontent.com/aaronmck/FlashFry/master/test_data/quickstart_data.tar.gz
tar xf quickstart_data.tar.gz

# testing
# database creation 
mkdir tmp
java -Xmx4g -jar /home/hx37930/software/FlashFry/FlashFry-assembly-1.15.jar \
index \
--tmpLocation ./tmp \
--database chr22_cas9ngg_database \
--reference chr22.fa.gz \
--enzyme spcas9ngg

# discover candidate targets and their potential off-target effets
java -Xmx4g -jar /home/hx37930/software/FlashFry/FlashFry-assembly-1.15.jar \
discover \
--database chr22_cas9ngg_database \
--fasta EMX1_GAGTCCGAGCAGAAGAAGAAGGG.fasta \
--output EMX1.output

# score the discovered sites
java -Xmx4g -jar /home/hx37930/software/FlashFry/FlashFry-assembly-1.15.jar \
score \
--input EMX1.output \
--output EMX1.output.scored \
--scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
--database chr22_cas9ngg_database

############################################################
# step 1: prepare target sequence files
############################################################
# assign input and output
dataDir="/scratch/hx37930/project/MTAG/06.crispr/data"
reference="/scratch/hx37930/project/MTAG/06.crispr/reference/hg38.fa"
database="/scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/hg38_cas9ngg_database"
outDir="/scratch/hx37930/project/MTAG/06.crispr/results/gRNAs"
inFile="${dataDir}/fine_mappedSNP.all.hg19Tohg38.map.txt"
output_prefix="targetRegion200bp_finemPIP05"

cd ${dataDir}

# create a bed file to extract a 200 bp region centering each fine-mapped SNPs
# create a bed file
awk 'BEGIN{FS=OFS="\t"}NR>1{print "chr"$1":"strtonum($2)-100"-"strtonum($2)+100}' ${inFile} > ${output_prefix}.bed
# extract target regions
ml SAMtools/1.17-GCC-12.2.0
samtools faidx $reference -r ${output_prefix}.bed -o  ${output_prefix}.fasta

############################################################
# step 2: database creation
############################################################
java -Xmx4g -jar /home/hx37930/software/FlashFry/FlashFry-assembly-1.15.jar \
        index \
        --tmpLocation ${outDir}/tmp \
        --database ${database} \
        --reference $reference \
        --enzyme spcas9ngg

############################################################
# step 3: discover candidate targets and their potential off-target effets
############################################################
java -Xmx4g -jar /home/hx37930/software/FlashFry/FlashFry-assembly-1.15.jar \
        discover \
        --database ${database} \
        --fasta ${dataDir}/${output_prefix}.fasta \
        --output ${outDir}/${output_prefix}.gRNAs.txt

############################################################
# step 4: score the discovered sites
############################################################
java -Xmx4g -jar /home/hx37930/software/FlashFry/FlashFry-assembly-1.15.jar \
	score \
	--input ${outDir}/${output_prefix}.gRNAs.txt \
	--output ${outDir}/${output_prefix}.gRNAs.scored.txt \
	--scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
	--database ${database}

############################################################
# step 5: combine candidate gRNAs with target finemapped SNPs using target regions
############################################################
##Usage: Rscript summary_gRNAs.r <length of target region /2> <scored file> <position file> <filter file> <output prefix> <tissue info: yes|no> <tissueCategory file if yes>
# all
Rscript summary_gRNAs.r 100 /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.all.gRNAs.scored.txt /scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.all.hg19Tohg38.map.txt /scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.all.removeTwoRegions.hg19Tohg38.map.txt /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.all.removeTwoRegions no
# 0.1 <= PIP < 0.5
Rscript summary_gRNAs.r 100 /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.all.gRNAs.scored.txt /scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.all.hg19Tohg38.map.txt /scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.PIP01-05.removeTwoRegions.sort.txt /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.PIP01-05.removeTwoRegions no
# PIP > 0.5
Rscript summary_gRNAs.r 100 /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.all.gRNAs.scored.txt /scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.all.hg19Tohg38.map.txt /scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.PIP05.removeTwoRegions.sort.txt /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.PIP05.removeTwoRegions yes /scratch/hx37930/project/MTAG/06.crispr/data/fine-mappedSNP.PIP05_tissueCategory.txt

############################################################
# step 5 (old version): combine candidate gRNAs with target finemapped SNPs using target regions
############################################################
# add target regions into mapped file
awk 'BEGIN{FS=OFS="\t"}{print $0,$1":"strtonum($2)-100"-"strtonum($2)+100}' ${inFile} > ${output_prefix}.hg19Tohg38.map.contigs.txt
awk -F '\t' 'NR==FNR{a=$5;b[a]=$1"\t"$2"\t"$3"\t"$4;next}{OFS="\t";c=$1;d[c]=$0;if(b[c]){print b[c],d[c]}}' ${output_prefix}.hg19Tohg38.map.contigs.txt ${outDir}/${output_prefix}.gRNAs.scored.txt > ${outDir}/${output_prefix}.gRNAs.scored.rsID.txt
awk -F '\t' 'NR==FNR{a=$5;b[a]=$1"\t"$2"\t"$3"\t"$4;next}{OFS="\t";c=$1;d[c]=$0;if(b[c]){print b[c],d[c]}}' ${output_prefix}.hg19Tohg38.map.contigs.txt ${outDir}/${output_prefix}.gRNAs.txt > ${outDir}/${output_prefix}.gRNAs.rsID.txt

# extract gRNAs for variants with PIP >= 0.5
head -n 1 ${outDir}/${output_prefix}.gRNAs.scored.rsID.txt > ${outDir}/fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.txt
awk -F '\t' 'NR==FNR{a=$1;b[a]=$0;next}{OFS="\t";c=$4;d[c]=$0;if(b[c]){print d[c]}}' fine_mappedSNP.PIP05.removeTwoRegions.sort.txt ${outDir}/${output_prefix}.gRNAs.scored.rsID.txt >> ${outDir}/fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.txt
head -n 1 ${outDir}/${output_prefix}.gRNAs.rsID.txt > ${outDir}/fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.rsID.txt
awk -F '\t' 'NR==FNR{a=$1;b[a]=$0;next}{OFS="\t";c=$4;d[c]=$0;if(b[c]){print d[c]}}' fine_mappedSNP.PIP05.removeTwoRegions.sort.txt ${outDir}/${output_prefix}.gRNAs.rsID.txt >> ${outDir}/fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.rsID.txt

# extract gRNAs for variants with 0.1 <= PIP < 0.5
head -n 1 ${outDir}/${output_prefix}.gRNAs.scored.rsID.txt > ${outDir}/fine_mappedSNP.PIP01-05.removeTwoRegions.gRNAs.scored.rsID.txt
awk -F '\t' 'NR==FNR{a=$1;b[a]=$0;next}{OFS="\t";c=$4;d[c]=$0;if(b[c]){print d[c]}}' fine_mappedSNP.PIP01-05.removeTwoRegions.sort.txt ${outDir}/${output_prefix}.gRNAs.scored.rsID.txt >> ${outDir}/fine_mappedSNP.PIP01-05.removeTwoRegions.gRNAs.scored.rsID.txt
head -n 1 ${outDir}/${output_prefix}.gRNAs.rsID.txt > ${outDir}/fine_mappedSNP.PIP01-05.removeTwoRegions.gRNAs.rsID.txt
awk -F '\t' 'NR==FNR{a=$1;b[a]=$0;next}{OFS="\t";c=$4;d[c]=$0;if(b[c]){print d[c]}}' fine_mappedSNP.PIP01-05.removeTwoRegions.sort.txt ${outDir}/${output_prefix}.gRNAs.rsID.txt >> ${outDir}/fine_mappedSNP.PIP01-05.removeTwoRegions.gRNAs.rsID.txt

# There are no gRNA results for some variants, extract those variants to design gRNAs using web server CHOPCHOP.
# check variants with PIP >= 0.5 but without gRNA results
cat fine_mappedSNP.PIP05.removeTwoRegions.sort.txt > ${outDir}/check_PIP05.scored.txt
awk 'NR>1{print $4}' ${outDir}/fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.txt >> ${outDir}/check_PIP05.scored.txt
sort ${outDir}/check_PIP05.scored.txt |uniq -c | awk '{if($1<2)print $2}' > ${outDir}/nogRNAs.tmp.txt
head -n 1 fine_mappedSNP.all.hg19Tohg38.map.contigs.txt > ${outDir}/nogRNAs.PIP05.removeTwoRegions.txt
awk -F '\t' 'NR==FNR{a=$4;b[a]=$0;next}{OFS="\t";c=$1;d[c]=$0;if(b[c]){print b[c]}}' fine_mappedSNP.all.hg19Tohg38.map.contigs.txt  ${outDir}/nogRNAs.tmp.txt >> ${outDir}/nogRNAs.PIP05.removeTwoRegions.txt
rm ${outDir}/check_PIP05.scored.txt ${outDir}/nogRNAs.tmp.txt
# check variants with 0.1 <= PIP < 0.5 but without gRNA results
cat fine_mappedSNP.PIP01-05.removeTwoRegions.sort.txt > ${outDir}/check_PIP01-05.scored.txt
awk 'NR>1{print $4}' ${outDir}/fine_mappedSNP.PIP01-05.removeTwoRegions.gRNAs.scored.rsID.txt >> ${outDir}/check_PIP01-05.scored.txt
sort ${outDir}/check_PIP01-05.scored.txt |uniq -c | awk '{if($1<2)print $2}' > ${outDir}/nogRNAs.tmp.txt
head -n 1 fine_mappedSNP.all.hg19Tohg38.map.contigs.txt > ${outDir}/nogRNAs.PIP01-05.removeTwoRegions.txt
awk -F '\t' 'NR==FNR{a=$4;b[a]=$0;next}{OFS="\t";c=$1;d[c]=$0;if(b[c]){print b[c]}}' fine_mappedSNP.all.hg19Tohg38.map.contigs.txt  ${outDir}/nogRNAs.tmp.txt >> ${outDir}/nogRNAs.PIP01-05.removeTwoRegions.txt
rm ${outDir}/check_PIP01-05.scored.txt ${outDir}/nogRNAs.tmp.txt


# merge variants with function annotations
cp fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.rsID.txt fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.rsID.tissue.txt
sed -i 's/#/\t/g' fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.rsID.tissue.txt
cp fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.txt fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.tissue.txt
sed -i 's/#/\t/g' fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.tissue.txt
cp nogRNAs.PIP05.removeTwoRegions.txt nogRNAs.PIP05.removeTwoRegions.tissue.txt
sed -i 's/#/\t/g' nogRNAs.PIP05.removeTwoRegions.tissue.txt
awk -F '\t' 'NR==FNR{a=$1;b[a]=$2;next}{OFS="\t";c=$4;d[c]=$0;if(b[c]){print d[c],b[c]}else {print d[c]"\tUnknown"}}' fine-mappedSNP_tissueCategory.txt  fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.rsID.tissue.txt > fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.rsID.tissue.txt2
mv fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.rsID.tissue.txt2 fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.rsID.tissue.txt
awk -F '\t' 'NR==FNR{a=$1;b[a]=$2;next}{OFS="\t";c=$4;d[c]=$0;if(b[c]){print d[c],b[c]}else {print d[c]"\tUnknown"}}' fine-mappedSNP_tissueCategory.txt fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.tissue.txt > fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.tissue.txt2
mv fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.tissue.txt2 fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.tissue.txt
awk -F '\t' 'NR==FNR{a=$1;b[a]=$2;next}{OFS="\t";c=$4;d[c]=$0;if(b[c]){print d[c],b[c]}else {print d[c]"\tUnknown"}}' fine-mappedSNP_tissueCategory.txt nogRNAs.PIP05.removeTwoRegions.tissue.txt > nogRNAs.PIP05.removeTwoRegions.tissue.txt2
mv nogRNAs.PIP05.removeTwoRegions.tissue.txt2 nogRNAs.PIP05.removeTwoRegions.tissue.txt

# extract top 3 gRNAs for each target variants
# sort file by chr, pos, Hsu2013 score
inFile="fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.tissue.txt"
outFile="fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.scored.rsID.tissue.top3.txt"
head -n 1 ${inFile} > ${outFile}
#awk 'NR>1{print $4}'  ${inFile} |sort |uniq > id.txt
awk 'NR>1{print $9}'  ${inFile} |sort |uniq > id.txt
for i in `cat id.txt`
do
	#grep $i ${inFile} |sort -r -k 14 |head -n 3 >> ${outFile}
	grep $i ${inFile} |sort -r -k 18 |head -n 3 >> ${outFile}
done
#awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5}' ${outFile} | uniq -c |awk '{if($1<3)print $2"\t"$3"\t"$4"\t"$5"\t"$6}' > lessthan3gRNAs.PIP01-05.removeTwoRegions.txt
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$29}' ${outFile} | uniq -c |awk '{if($1<3)print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}'  > lessthan3gRNAs.PIP05.removeTwoRegions.tissue.txt
