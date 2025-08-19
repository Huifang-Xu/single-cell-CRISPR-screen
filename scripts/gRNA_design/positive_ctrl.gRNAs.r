#!Rscript

library(dplyr)
library(tidyverse)
library(readr)
library(data.table)
library(biomaRt)

# read in FADS knockout experimental data (three replicates)
rep1 <- fread("/scratch/hx37930/project/MTAG/06.crispr/data/positive_control_FADS_08082024/22102R-01-01_S22_.txt",header=T,sep="\t")

rep2 <- fread("/scratch/hx37930/project/MTAG/06.crispr/data/positive_control_FADS_08082024/22102R-01-02_S23_.txt",header=T,sep="\t")

rep3 <- fread("/scratch/hx37930/project/MTAG/06.crispr/data/positive_control_FADS_08082024/22102R-01-03_S24_.txt",header=T,sep="\t")

effect.genes <- fread("/scratch/hx37930/project/MTAG/06.crispr/data/positive_control_FADS_08082024/gene_effects_ACH-000739.csv",header=T,sep=",")

# merge three replicates
three.reps <- full_join(rep1, rep2, by="Transcript_Id")
three.reps <- full_join(three.reps, rep3, by="Transcript_Id")

# calculate average gene expression level
three.reps$mean.TPM <- rowMeans(three.reps[,2:4])

# order by highest average gene expression level
three.reps.order <- three.reps[order(three.reps$mean.TPM,decreasing=TRUE),]

# Retrieve gene names from transcript IDs
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
transcript_ids <- three.reps.order$Transcript_Id  
gene_info <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name","chromosome_name", "start_position","end_position", "strand"),
                   filters = "ensembl_transcript_id",
                   values = transcript_ids,
                   mart = ensembl)
three.reps.gene_info <- left_join(three.reps.order,gene_info,by=c("Transcript_Id"="ensembl_transcript_id"))

# add gene effects
three.reps.geneEffects <- left_join(three.reps.gene_info,effect.genes,by=c("external_gene_name"="gene"))
# Retrieve non-essential genes
three.reps.nonessen <- three.reps.geneEffects[three.reps.geneEffects$z_score > 0, ]

# 20 non-esstential genes
trascript.names <- c("ENST00000484695","ENST00000619601","ENST00000489769","ENST00000551679","ENST00000550343",
		     "ENST00000273784","ENST00000583217","ENST00000474253","ENST00000565355","ENST00000520028",
		     "ENST00000464611","ENST00000596179","ENST00000505309","ENST00000480603","ENST00000474281",
		     "ENST00000435330","ENST00000526862","ENST00000546656","ENST00000487713","ENST00000468955")
gene.names <- c("FGG","GAPDH","SERPINA1","HNRNPA1","NACA",
		"AHSG","EIF4A1","TPI1","ALDOA","HINT1",
		"ACTB","C3","SELENOP","PPIA","PGK1",
		"RPL7","HSPA8","KRT18","SAT1","IGFBP1")

		     
# Retrieve top 20 non-essential genes highly expressed in HepG2
three.reps.nonessen.top20 <- subset(three.reps.nonessen,Transcript_Id %in% trascript.names)

# add FADS1 and FADS2 into positive controls
three.reps.nonessen.top20.fads <- rbind(three.reps.nonessen.top20,three.reps.geneEffects[three.reps.geneEffects$Transcript_Id=="ENST00000350997" | three.reps.geneEffects$Transcript_Id=="ENST00000278840",])
 
# add contig info: +- 200bp of TSS, or upstram 400bp of TSS
# Function to generate the 5th column output
add_contig <- function(chr, start, end, strand,length) {
  if (strand == -1) {
    paste0("chr",chr, ":", (end - length), "-", (end + length))
  } else if (strand == 1) {
    paste0("chr",chr, ":", (start - length), "-", (start + length))
  } else {
    NA  # In case there's an unexpected value in the strand column
  }
}

# Apply the function to each row to create a new column "contig": 200bp, 400bp, 800bp
three.reps.nonessen.top20.fads$contig.100bp <- mapply(add_contig, three.reps.nonessen.top20.fads$chromosome_name, three.reps.nonessen.top20.fads$start_position, three.reps.nonessen.top20.fads$end_position,three.reps.nonessen.top20.fads$strand,100)
three.reps.nonessen.top20.fads$contig.200bp <- mapply(add_contig, three.reps.nonessen.top20.fads$chromosome_name, three.reps.nonessen.top20.fads$start_position, three.reps.nonessen.top20.fads$end_position,three.reps.nonessen.top20.fads$strand,200)
three.reps.nonessen.top20.fads$contig.400bp <- mapply(add_contig, three.reps.nonessen.top20.fads$chromosome_name, three.reps.nonessen.top20.fads$start_position, three.reps.nonessen.top20.fads$end_position,three.reps.nonessen.top20.fads$strand,400)

# output table
write.table(three.reps.nonessen.top20.fads,file="/scratch/hx37930/project/MTAG/06.crispr/data/positive_control_FADS_08082024/positive_control_22genes.txt",quote=F,col.names=T,row.names=F,sep="\t")
