#!Rscript

##Usage: Rscript summary_gRNAs.r <length of target region/2> <scored file> <position file> <filter file> <output prefix> <tissue info: yes|no> <tissueCategory file>
# Example: Rscript summary_gRNAs.r 100 /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.all.gRNAs.scored.txt /scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.all.hg19Tohg38.map.txt /scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.PIP05.removeTwoRegions.sort.txt /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.PIP05.removeTwoRegions yes /scratch/hx37930/project/MTAG/06.crispr/data/fine-mappedSNP.PIP05_tissueCategory.txt

library(dplyr)
library(tidyverse)
library(readr)
library(data.table)

# accept parameters
Args <- commandArgs()

length_target <- Args[6]

# read in data
#score.all.file <- "/scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.all.gRNAs.scored.txt"
score.all.file <- Args[7]
score.all <- fread(score.all.file,header=T, sep="\t")

#position.file <- "/scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.all.hg19Tohg38.map.txt"
position.file <- Args[8]
position.all <- fread(position.file,header=T, sep="\t")

#filter.file <- "/scratch/hx37930/project/MTAG/06.crispr/data/fine_mappedSNP.PIP05.removeTwoRegions.sort.txt"
filter.file <- Args[9]
filter.set <- fread(filter.file,header=T, sep="\t")

#output.prefix <- "/scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/fine_mappedSNP.PIP05.removeTwoRegions"
output.prefix <- Args[10] 

#tissueInfo <- "no"
tissueInfo <- Args[11]
if (tissueInfo == "yes") {
#	tissue.file <- "/scratch/hx37930/project/MTAG/06.crispr/data/fine-mappedSNP.PIP05_tissueCategory.txt"
	tissue.file <- Args[12]
	tissue <- fread(tissue.file,header=T, sep="\t")
}

# add contig info to position file: 400bp (+- 200), 200bp (+- 100)
position.all$contig <- paste(position.all$hg38chr,":",as.numeric(position.all$hg38start)-as.numeric(length_target),"-",as.numeric(position.all$hg38start)+as.numeric(length_target),sep="")

# combine score file with position file
score.comb <- left_join(position.all,score.all,by="contig") 

# extract variants in the filter set
score.filter <- subset(score.comb, `SNP:hg19chr:hg19BP:A1:A2` %in% filter.set$`SNP:hg19chr:hg19BP:A1:A2`)

# filter gRNAs based on "Hsu2013" score >= 50
score.filter.offtarget <- score.filter[score.filter$Hsu2013 >= 50,]

# extract variants with > 3 gRNAs and Hsu2013 score > 50
gRNAs.count <- as.data.frame(table(score.filter.offtarget$contig))
gt3grnas <- gRNAs.count[gRNAs.count$Freq >=3,"Var1"]
score.filter.offtarget.3gRNAs <- subset(score.filter.offtarget, contig %in% gt3grnas)
score.filter.offtarget.order <- score.filter.offtarget.3gRNAs[with(score.filter.offtarget.3gRNAs,order(contig,-Hsu2013)),]
score.filter.offtarget.top3 <- score.filter.offtarget.order %>% group_by(contig) %>% slice_head(n = 3) %>% ungroup()

# For some variants, less than 3 gRNAs were produced and needed to be re-design. Create a list to store these variants.
redesign.snp <- unique(subset(score.filter[,c("hg38chr","hg38start","hg38end","SNP:hg19chr:hg19BP:A1:A2")],!(`SNP:hg19chr:hg19BP:A1:A2` %in% score.filter.offtarget.top3$`SNP:hg19chr:hg19BP:A1:A2`)))

# if tissue information is provided, then combine score file with tissue information
if (tissueInfo == "yes") {
	# split column "SNP:hg19chr:hg19BP:A1:A2" to get rsID
	score.filter.offtarget.top3.split <- score.filter.offtarget.top3
	score.filter.offtarget.top3.split$uniqID <- score.filter.offtarget.top3.split$`SNP:hg19chr:hg19BP:A1:A2`
	score.filter.offtarget.top3.split <- separate_wider_delim(score.filter.offtarget.top3.split, cols = uniqID, delim = "#", names = c("SNP", "hg19chr","hg19BP","A1","A2"))
	# combine tissue info with score file
	score.filter.offtarget.top3.tissue <- left_join(score.filter.offtarget.top3.split,tissue,by="SNP")
	#subset required columns
	colnames.all <- c(names(score.filter.offtarget.top3),names(tissue)[2])
	score.tissue.result <- score.filter.offtarget.top3.tissue[,colnames.all]
	score.tissue.result.file <- paste(output.prefix,".gRNAs.scored.rsID.top3.tissue.txt",sep="")
	write.table(score.tissue.result,file=score.tissue.result.file,quote=F,col.names=T,row.names=F,sep="\t")
}

# output tables
redesign.snp.file <- paste(output.prefix,"gRNAslt3.txt",sep="")
write.table(redesign.snp,file=redesign.snp.file,quote=F,col.names=T,row.names=F,sep="\t")
score.filter.offtarget.top3.file <- paste(output.prefix,".gRNAs.scored.rsID.top3.txt",sep="")
write.table(score.filter.offtarget.top3,file=score.filter.offtarget.top3.file,quote=F,col.names=T,row.names=F,sep="\t")
