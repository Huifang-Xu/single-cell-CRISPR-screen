#!bin/Rscript
# R version 4.3.1 R/4.3.1-foss-2022a

# load libraries
library(optparse)

# define parameters
option_list <- list(
	make_option(c("-o", "--outputDir"), type = 'character', default = getwd(), help = "output directory [default = %default]", metavar = "character"),
	make_option(c("-g", "--geneInfo"), type = 'character', default=NULL, help = "gene information (e.g. all_genes.csv)", metavar = "CHARACTER"),
	make_option(c("-r", "--guideInfo"), type = 'character', default=NULL, help = "guide information (e.g. all_guides.clean.csv)", metavar = "CHARACTER"),
	make_option(c("-G", "--geneMatrix"), type = 'character', default=NULL, help = "gene count matrix (e.g. DGE.mtx)", metavar = "CHARACTER"),
	make_option(c("-R", "--guideMatrix"), type = 'character', default=NULL, help = "guide RNA count matrix (e.g. count_matrix.mtx)", metavar = "CHARACTER"),
	make_option(c("-c", "--guideChr"), type = 'character', default=NULL, help = "guide information adding chr, position (e.g. guideList.chr.uniq.txt)", metavar = "CHARACTER"),
	make_option(c("-i", "--moi"), type = 'character', default = "high", help = "multiplicity of infection (MOI) [default = %default]", metavar = "CHARACTER"),
	make_option(c("-m", "--mitoProp"), type = 'double', default = 0.075, help = "cutoff to remove cells with high fraction of mitochondrial genes [default = %default]", metavar = "NUMBER"),
	make_option(c("-t", "--type"), type = 'character', default = "cis", help = "run mode: cis, trans, or both [default = %default]", metavar = "CHARACTER")
	)
# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
Args <- parse_args(opt_parser)
print(Args)
#Args <- commandArgs(trailingOnly=TRUE)
# Check if arguments are provided
if (length(Args) == 0) {
  stop("No arguments provided")
}

directory_to_write <- Args$outputDir 
all_genes_fp <- Args$geneInfo
all_grnas_fp <- Args$guideInfo
gene_mat_fp <- Args$geneMatrix
grna_mat_fp <- Args$guideMatrix
grna_target_data_frame <- Args$guideInfo
guide_chr <- Args$guideChr
moi <- Args$moi
p_mito_threshold <- Args$mitoProp
type <- Args$type

#### example inputs and outputs
# assign inputs and outputs
#directory_to_write <- "/scratch/../sceptre/yk3_19117_35"
# genes
#all_genes_fp <- "/scratch/../combined_parents/yk3/DGE_filtered/all_genes.csv"
# gRNAs
#all_grnas_fp <- "/scratch/../combined_crispr/yk3/guide_RNAs_filtered/all_guides.clean.csv"
# gene matrix
#gene_mat_fp <- "/scratch/../combined_parents/yk3/DGE_filtered/DGE.mtx"
# gRNAs matrix
#grna_mat_fp <- "/scratch/../combined_crispr/yk3/guide_RNAs_filtered/count_matrix.mtx"
# grna-target info
#grna_target_data_frame <- "/scratch/../combined_crispr/yk3/guide_RNAs_filtered/all_guides.clean.csv"
#guide_chr <- "/scratch/../yk3_19117_35/guideList.chr.uniq.txt"
#moi <- "high" # moi: high or low
#p_mito_threshold <- 0.2 # mitochrondria genes cutoff

# load libraries
library(sceptre)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)


# Convenience functions for plotting
SaveFigure <- function(plots, name, type = "png", width, height, res){
        if(type == "png") {
                png(paste0(directory_to_write, "/", name, ".", type), type=c("cairo", "cairo-png", "Xlib", "quartz"), width = width, height = height, units = "in", res = 300)
        } else {
                pdf(paste0(directory_to_write, "/", name, ".", type), width = width, height = height)
        }
        print(plots)
        dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(directory_to_write, "/", name, ".rds"))
}

ReadObject <- function(name){
  readRDS(paste0(directory_to_write, "/", name, ".rds"))
}

# construct positive control pairs: 1) set grna_id and grna_target, 2) change the grna_target of "ntc" to "non-targeting", change gene name containing "_1/_2/_3" to original gene name, 3) get ensmbl IDs
grna_target_data_frame <- fread(grna_target_data_frame)
grna_target_data_frame <- grna_target_data_frame[,c(2,2,3)]
colnames(grna_target_data_frame) <- c("grna_id","grna_target","genome")
grna_target_data_frame$grna_target[grepl("ntc",grna_target_data_frame$grna_target)] <- "non-targeting"

positive_control_grnas <- subset(grna_target_data_frame,! (grepl("Liver",grna_target) | (grepl("non-targeting",grna_target))))
positive_control_grnas$grna_target <- as.data.frame(apply(as.data.frame(positive_control_grnas$grna_target), 2, function(x){gsub(pattern = "_1", replacement = "",x) }))
positive_control_grnas$grna_target <- as.data.frame(apply(as.data.frame(positive_control_grnas$grna_target), 2, function(x){gsub(pattern = "_2", replacement = "",x) }))
positive_control_grnas$grna_target <- as.data.frame(apply(as.data.frame(positive_control_grnas$grna_target), 2, function(x){gsub(pattern = "_3", replacement = "",x) }))
#names(positive_control_grnas) <- c("grna_id","grna_target","genome")

# change response_id to ensembl "ENSG" id
gene_names <- as.vector(unique(positive_control_grnas$grna_target))
gene_names_ids <- as.data.frame(mapIds(org.Hs.eg.db,
    keys = gene_names,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"))
colnames(gene_names_ids) <- "response_id"
gene_names_ids$gene_name <- rownames(gene_names_ids)

positive_control_pairs_merge <- left_join(positive_control_grnas,gene_names_ids,by=c("grna_target"="gene_name"))
positive_control_pairs <- positive_control_pairs_merge[,c("grna_id","response_id")]
colnames(positive_control_pairs) <- c("grna_target","response_id")
#positive_control_pairs$grna_group <- str_split_fixed(positive_control_pairs$grna_target,"_",2)[,1]


# negative control pairs
negative_control_pairs <- subset(grna_target_data_frame[,c("grna_id","grna_target","grna_id")], grepl("non-targeting",grna_target))
colnames(negative_control_pairs) <- c("grna_id","grna_target","grna_group")


grna_target_data_frame <- as.data.frame(grna_target_data_frame)
guide_chr <- fread(guide_chr)
grna_target_data_frame <- left_join(grna_target_data_frame,guide_chr[,c("Guide_Name","chr","start","end")],by=c("grna_id"="Guide_Name"))

# add causal variants and targeting region of each gRNA
guideList_causalSNPs <- fread("/scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/guideList.contig.uniq.txt",header=T)

# create a sceptre_object by importing data generated from ParseBio
sceptre_object <- import_data_from_parse(
  gene_mat_fp = gene_mat_fp,
  grna_mat_fp = grna_mat_fp,
  all_genes_fp = all_genes_fp,
  all_grnas_fp = all_grnas_fp,
  moi = moi,
  grna_target_data_frame = grna_target_data_frame)
# check slots: slotNames(sceptre_object)
sceptre_object@negative_control_pairs <- negative_control_pairs

SaveObject(sceptre_object,"sceptre_object.beforeQC")

# response matrix: return # of features, # of cells, backing file path
response_matrix <- get_response_matrix(sceptre_object)
dim(response_matrix)
# [1]  62710 246793

response_ids <- rownames(response_matrix)

# gRNA matrix: return # of features, # of cells, backing file path
grna_matrix <- get_grna_matrix(sceptre_object)
dim(grna_matrix)
# [1]   1184 246793

cell_covariates <- get_cell_covariates(sceptre_object)
# add extra covariates: batch

write.table(cell_covariates,paste0(directory_to_write,"/cell_covariates.txt"),quote=T,row.names=F,col.names=T)

# construct_cis_pairs(): returns the set of target-response pairs for which the target and response are located on the same chromosome and in close physical proximity to one another (default 500Kb).
discovery_pairs_cis <- construct_cis_pairs(sceptre_object,positive_control_pairs = positive_control_pairs)
# distance_threshold = 5e6, # 5Mb

# exclude all pairs containing a positive control gRNA target from the trans pairs
discovery_pairs_trans <- construct_trans_pairs(sceptre_object = sceptre_object, positive_control_pairs = positive_control_pairs,pairs_to_exclude = "pairs_containing_pc_targets")

# cis: side="left", grna_integration_strategy = "union", control="complement set"
#trans: side="both", grna_integration_strategy = "union", control="complement set"

directory_to_write_old <- directory_to_write
######################################################
############### run cis analysis ###################
######################################################
if (type == "cis" | type == "both") {
	if (type == "both") {
		# rename directory_to_write
		directory_to_write <- paste0(directory_to_write_old,"/cis")
	}

# set parameters
sceptre_object <- sceptre_object |>
  set_analysis_parameters(
    discovery_pairs = discovery_pairs_cis,
    positive_control_pairs = positive_control_pairs,
    side = "left",
    control_group = "complement")

#sceptre_object@positive_control_pairs <- positive_control_pairs

# assign gRNAs
sceptre_object <- sceptre_object |>  assign_grnas(parallel = TRUE)

# plot covariates before QC
plot.cis.covariates <- plot_covariates(sceptre_object, p_mito_threshold = p_mito_threshold)
SaveFigure(plot.cis.covariates, "dist_cis_covariates.beforeQC", width = 8, height = 6, res = 300)

# run QC
sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = p_mito_threshold)

SaveObject(sceptre_object,"sceptre_object.afterQC")

# plot covariates after QC
#plot.cis.covariates.qc <- plot_covariates(sceptre_object, p_mito_threshold = p_mito_threshold)
#SaveFigure(plot.cis.covariates.qc, "dist_cis_covariates.afterQC", width = 10, height = 6, res = 300)

# plot QC results
#plot.cis.qc <- plot(sceptre_object)
#SaveFigure(plot.cis.qc, "results_cis_qc",width = 6, height = 8,res=300)

# run calibration check
sceptre_object <- sceptre_object |> run_calibration_check(parallel = FALSE)

SaveObject(sceptre_object,"sceptre_object.calibration.cis")

# plot caliration results
plot.cis.calibration <- plot(sceptre_object)
SaveFigure(plot.cis.calibration,"results_cis_calibration", width = 10, height = 10, res =300)

# get calibration result
calibration_result.cis <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_calibration_check")

# run power check
sceptre_object <- sceptre_object |> run_power_check(parallel = TRUE)

# run discovery analysis
sceptre_object <- sceptre_object |> run_discovery_analysis(parallel = TRUE)

discovery_result.cis <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_discovery_analysis"
)

# get p-values of positive controls, negative controls and discovery pairs
power_result <- sceptre_object@power_result
power_result$significant <- NA
power_result$group <- "Positive control"
calibration_result.cis$group <- "Negative control"
discovery_result.cis$group <- "Discovery"
# combine results of positive, negative controls, and discovery
all_results <- rbind(power_result,calibration_result.cis,discovery_result.cis)

# add gene symbols based on ensembl IDs
ensembl_ids <- as.vector(unique(all_results$response_id))
# mapping IDs
gene_symbol <- mapIds(org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first")
gene_symbols_df <- data.frame(gene_symbol)
gene_symbols_df$response_id <- names(gene_symbol)

# update calibration_result
calibration_result.cis.updated <- left_join(calibration_result.cis,gene_symbols_df,by="response_id")
write.table(calibration_result.cis.updated,paste0(directory_to_write,"/calibration_result.cis.txt"),quote=F,row.names=F,col.names=T)

# update discovery_result: add gene symbol and targeting causal variant
discovery_result.cis.updated <- left_join(discovery_result.cis,gene_symbols_df,by="response_id")
discovery_result.cis.updated <- left_join(discovery_result.cis.updated,guideList_causalSNPs,by="grna_target")
write.table(discovery_result.cis.updated,paste0(directory_to_write,"/discovery_result.cis.txt"),quote=F,row.names=F,col.names=T)

# update all_results
all_results <- left_join(all_results,gene_symbols_df,by="response_id")
all_results <- left_join(all_results,guideList_causalSNPs,by="grna_target")
write.table(all_results,paste0(directory_to_write,"/all_result.cis.txt"),quote=F,row.names=F,col.names=T,sep="\t")

# write outputs
write_outputs_to_directory(
  sceptre_object = sceptre_object,
  directory = directory_to_write)

# construct obeserved- and expected- p-values
qq_ptc <- data.frame(obs=-log10(sort(power_result$p_value)))
qq_ptc$exp <- -log10(seq(length(qq_ptc$obs))/(length(qq_ptc$obs)+1))
qq_ptc$group <- "Positive control"
qq_ntc <- data.frame(obs=-log10(sort(calibration_result.cis$p_value)))
qq_ntc$exp <- -log10(seq(length(qq_ntc$obs))/(length(qq_ntc$obs)+1))
qq_ntc$group <- "Negative control"
qq_discovery <- data.frame(obs=-log10(sort(discovery_result.cis$p_value)))
qq_discovery$exp <- -log10(seq(length(qq_discovery$obs))/(length(qq_discovery$obs)+1))
qq_discovery$group <- "Discovery"
qq <- rbind(qq_ptc,qq_ntc,qq_discovery)

#plot QQ plot
library(RColorBrewer)
colors <- brewer.pal(3,"Dark2")
plot_qq <- ggplot(qq,aes(exp,obs)) + geom_point(aes(colour=group))+ geom_abline(intercept=0,slope=1,colour="lightblue")+
scale_x_continuous(expression(paste("Expected ",-log[10],"P")),expand=c(0,0.1))+scale_y_continuous(expression(paste("Observed ",-log[10],"P")),expand=c(0,0.1))+theme(panel.grid=element_blank(),panel.background=element_blank(),axis.line=element_line(),axis.text=element_text(colour="black"),legend.position=c(0.55,0.88),legend.title=element_blank())+scale_color_manual(values = colors)
SaveFigure(plot_qq, "plot_qq_cis",width = 5, height = 5,res=300)
}

######################################################
############### run trans analysis ###################
######################################################
if (type == "trans" | type == "both") {
	if (type == "both") {
                # rename directory_to_write
                directory_to_write <- paste0(directory_to_write_old,"/trans")
        }

sceptre_object <- sceptre_object |>
  set_analysis_parameters(
    discovery_pairs = discovery_pairs_trans,
    positive_control_pairs = positive_control_pairs,
    side = "both",
    control_group = "complement") 

# assign gRNAs
sceptre_object <- sceptre_object |>  assign_grnas(parallel = TRUE) 

# plot covariates before QC
plot.trans.covariates <- plot_covariates(sceptre_object, p_mito_threshold = p_mito_threshold)
SaveFigure(plot.trans.covariates, "dist_trans_covariates.beforeQC", width = 10, height = 6, res = 300)

# run QC 
sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = p_mito_threshold)

SaveObject(sceptre_object,"sceptre_object.afterQC")

# plot covariates after QC
#plot.trans.covariates.qc <- plot_covariates(sceptre_object, p_mito_threshold = p_mito_threshold)
#SaveFigure(plot.trans.covariates.qc, "dist_trans_covariates.afterQC", width = 10, height = 6, res = 300)

# plot QC results
#plot.trans.qc <- plot(sceptre_object)
#SaveFigure(plot.trans.qc, "results_trans_qc",width = 6, height = 8,res=300)

# run calibration check
sceptre_object <- sceptre_object |> run_calibration_check(parallel = TRUE)
  
SaveObject(sceptre_object,"sceptre_object.calibration.trans")

# plot caliration results
#plot.trans.calibration <- plot(sceptre_object)
#SaveFigure(plot.trans.calibration,"results_trans_calibration", width = 10, height = 10, res =300)

# get calibration result
calibration_result.trans <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_calibration_check")

# run power check
sceptre_object <- sceptre_object |> run_power_check(parallel = TRUE) 

# run discovery analysis
sceptre_object <- sceptre_object |> run_discovery_analysis(parallel = TRUE)

discovery_result.trans <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_discovery_analysis"
)

write_outputs_to_directory(
  sceptre_object = sceptre_object,
  directory = directory_to_write)

# get p-values of positive controls
power_result.trans <- sceptre_object@power_result
power_result.trans$significant <- NA
power_result.trans$group <- "Positive control"
calibration_result.trans$group <- "Negative control"
discovery_result.trans$group <- "Discovery"
# combine results of positive, negative controls, and discovery
all_results.trans <- rbind(power_result.trans,calibration_result.trans,discovery_result.trans)

# add gene symbols based on ensembl IDs
ensembl_ids <- as.vector(unique(all_results.trans$response_id))
# mapping IDs
gene_symbol <- mapIds(org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first")
gene_symbols_df <- data.frame(gene_symbol)
gene_symbols_df$response_id <- names(gene_symbol)

# update calibration_result
calibration_result.trans.updated<- left_join(calibration_result.trans,gene_symbols_df,by="response_id")
write.table(calibration_result.trans.updated,paste0(directory_to_write,"/calibration_result.trans.txt"),quote=F,row.names=F,col.names=T)

# update discovery_result: add gene symbol and targeting causal variant
discovery_result.trans.updated <- left_join(discovery_result.trans,gene_symbols_df,by="response_id")
discovery_result.trans.updated <- left_join(discovery_result.trans.updated,guideList_causalSNPs,by="grna_target")
write.table(discovery_result.trans.updated,paste0(directory_to_write,"/discovery_result.trans.txt"),quote=F,row.names=F,col.names=T)

# update all_results.trans
all_results.trans <- left_join(all_results.trans,gene_symbols_df,by="response_id")
all_results.trans <- left_join(all_results.trans,guideList_causalSNPs,by="grna_target")
write.table(all_results.trans,paste0(directory_to_write,"/all_result.trans.txt"),quote=F,row.names=F,col.names=T,sep="\t")

#write_outputs_to_directory(
#  sceptre_object = sceptre_object,
#  directory = directory_to_write)

# construct obeserved- and expected- p-values
qq_ptc.trans <- data.frame(obs=-log10(sort(power_result.trans$p_value)))
qq_ptc.trans$exp <- -log10(seq(length(qq_ptc.trans$obs))/(length(qq_ptc.trans$obs)+1))
qq_ptc.trans$group <- "Positive control"
qq_ntc.trans <- data.frame(obs=-log10(sort(calibration_result.trans$p_value)))
qq_ntc.trans$exp <- -log10(seq(length(qq_ntc.trans$obs))/(length(qq_ntc.trans$obs)+1))
qq_ntc.trans$group <- "Negative control"
qq_discovery.trans <- data.frame(obs=-log10(sort(discovery_result.trans$p_value)))
qq_discovery.trans$exp <- -log10(seq(length(qq_discovery.trans$obs))/(length(qq_discovery.trans$obs)+1))
qq_discovery.trans$group <- "Discovery"
qq.trans <- rbind(qq_ptc.trans,qq_ntc.trans,qq_discovery.trans)

#plot QQ plot
library(RColorBrewer)
colors <- brewer.pal(3,"Dark2")
plot_qq.trans <- ggplot(qq.trans,aes(exp,obs)) + geom_point(aes(colour=group))+ geom_abline(intercept=0,slope=1,colour="lightblue")+
scale_x_continuous(expression(paste("Expected ",-log[10],"P")),expand=c(0,0.1))+scale_y_continuous(expression(paste("Observed ",-log[10],"P")),expand=c(0,0.1))+theme(panel.grid=element_blank(),panel.background=element_blank(),axis.line=element_line(),axis.text=element_text(colour="black"),legend.position=c(0.55,0.88),legend.title=element_blank())+scale_color_manual(values = colors)
SaveFigure(plot_qq.trans, "plot_qq_trans",width = 5, height = 5,res=300)
}
