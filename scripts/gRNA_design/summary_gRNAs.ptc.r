#!Rscript

##Usage: Rscript summary_gRNAs.r <scored file> <position file> <column # of contig> <output prefix>
# Example: Rscript summary_gRNAs.r /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/positive_control_FADS_08082024/positive_control_22genes.100bp.gRNAs.scored.txt /scratch/hx37930/project/MTAG/06.crispr/data/positive_control_FADS_08082024/positive_control_22genes.txt 15 /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/positive_control_FADS_08082024/positive_control_22genes.100bp

library(dplyr)
library(tidyverse)
library(readr)
library(data.table)

# accept parameters
Args <- commandArgs()

# read in data
score.all.file <- Args[6]
score.all <- fread(score.all.file,header=T, sep="\t")

position.file <- Args[7]
position.all <- fread(position.file,header=T, sep="\t")

column.contig <- Args[8]

output.prefix <- Args[9] 

# rename contig column
raw_name <- colnames(position.all)[column.contig]
names(position.all)[names(position.all)==raw_name] <- "contig"

# combine score file with position file
score.comb <- left_join(position.all,score.all,by="contig") 

# filter gRNAs based on "Hsu2013" score >= 50
score.filter.offtarget <- score.comb[score.comb$Hsu2013 >= 50,]

# extract variants with > 3 gRNAs and Hsu2013 score > 50
gRNAs.count <- as.data.frame(table(score.filter.offtarget$contig))
gt3grnas <- gRNAs.count[gRNAs.count$Freq >=3,"Var1"]
score.filter.offtarget.3gRNAs <- subset(score.filter.offtarget, contig %in% gt3grnas)
score.filter.offtarget.order <- score.filter.offtarget.3gRNAs[with(score.filter.offtarget.3gRNAs,order(contig,-Hsu2013)),]
score.filter.offtarget.top3 <- score.filter.offtarget.order %>% group_by(contig) %>% slice_head(n = 3) %>% ungroup()

# remove non-informative columns
keep.columns <- c("external_gene_name","Transcript_Id","mean.TPM","chromosome_name","start_position","end_position","mean",names(score.all))
score.filter.offtarget.top3.final <- score.filter.offtarget.top3[,keep.columns]

# For some variants, less than 3 gRNAs were produced and needed to be re-design. Create a list to store these variants.
redesign.snp <- unique(subset(position.all,!(contig %in% score.filter.offtarget.top3$contig)))

# output tables
if (nrow(redesign.snp) >=1 ) {
	redesign.snp.file <- paste(output.prefix,"gRNAslt3.txt",sep="")
	write.table(redesign.snp,file=redesign.snp.file,quote=F,col.names=T,row.names=F,sep="\t")
}

score.filter.offtarget.top3.file <- paste(output.prefix,".gRNAs.scored.rsID.top3.txt",sep="")
write.table(score.filter.offtarget.top3.final,file=score.filter.offtarget.top3.file,quote=F,col.names=T,row.names=F,sep="\t")
