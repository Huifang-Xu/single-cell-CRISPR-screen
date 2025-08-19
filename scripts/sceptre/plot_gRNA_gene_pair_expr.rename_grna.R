#!bin/Rscript

# ml R/4.3.1-foss-2022a

# load libraries
library(sceptre)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(readr)
library(tidyr)

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

# read in files
directory_to_write <- "/scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/merged_batch12/v3_rename_gRNA/cis"
sceptre_object <- ReadObject("sceptre_object.calibration.cis")
all_result <- fread(paste0(directory_to_write,"/all_result.cis.txt"))

positive_control_pairs <- all_result[all_result$group=="Positive control",]
negative_control_pairs <- all_result[all_result$group=="Negative control",]
discovery_pairs <- all_result[all_result$group=="Discovery",]

# plot positive controls
response_id_to_plot <- positive_control_pairs$response_id
for (i in 1:length(response_id_to_plot)) {
        grna_target_to_plot <- positive_control_pairs$grna_target[positive_control_pairs$response_id==response_id_to_plot[i]]
        gene_name <- grna_target_to_plot
        p_value <- sprintf(positive_control_pairs$p_value[positive_control_pairs$response_id==response_id_to_plot[i]],fmt = "%0.02e")
	logfc = sprintf(positive_control_pairs$log_2_fold_change[positive_control_pairs$response_id==response_id_to_plot[i]],fmt = "%.2f")
	plot.pairs <- plot_response_grna_target_pair(sceptre_object,response_id = response_id_to_plot[i],grna_target = as.character(grna_target_to_plot))
	plot.pairs <- annotate_figure(plot.pairs, top = text_grob(paste0(gene_name,"\n","p-value = ",p_value, "; log FC = ", logfc), face = "bold", size = 12))
        SaveFigure(plot.pairs,paste0("positive_pairs.",gene_name), width = 5, height = 5, res =300)
}

# plot discovery pairs
df_to_plot <- discovery_pairs[1:100,]
for (i in 1:nrow(df_to_plot)) {
	response_id_to_plot <- as.character(df_to_plot[i,"response_id"])
	snp_to_plot <- as.character(df_to_plot[i,"SNP"])
	grna_target_to_plot <- as.character(df_to_plot[i,'grna_target'])
	gene_name <- as.character(df_to_plot[i,"gene_symbol"])
	p_value <- sprintf(df_to_plot[i,"p_value"],fmt = "%0.02e")
	logfc = sprintf(df_to_plot[i,"log_2_fold_change"],fmt = "%.2f")
        plot.pairs <- plot_response_grna_target_pair(sceptre_object,response_id = response_id_to_plot,grna_target = as.character(grna_target_to_plot))
        plot.pairs <- annotate_figure(plot.pairs, top = text_grob(paste0(gene_name," : ",snp_to_plot, "\n","p-value = ", p_value, "; log FC = ", logfc), face = "bold", size = 12))
	SaveFigure(plot.pairs,paste0("discovery_pairs.",snp_to_plot,"_",gene_name), width = 5, height = 5, res =300)
}

# plot highlighting discovery pairs
cre_to_plot <- "target_rs28690720"
df_to_plot <- discovery_pairs[discovery_pairs$grna_target == cre_to_plot,]
for (i in 1:nrow(df_to_plot)) {
        response_id_to_plot <- as.character(df_to_plot[i,"response_id"])
        snp_to_plot <- as.character(df_to_plot[i,"SNP"])
        grna_target_to_plot <- as.character(df_to_plot[i,'grna_target'])
        gene_name <- as.character(df_to_plot[i,"gene_symbol"])
        p_value <- sprintf(df_to_plot[i,"p_value"],fmt = "%0.02e")
        logfc = sprintf(df_to_plot[i,"log_2_fold_change"],fmt = "%.2f")
        plot.pairs <- plot_response_grna_target_pair(sceptre_object,response_id = response_id_to_plot,grna_target = as.character(grna_target_to_plot))
        plot.pairs <- annotate_figure(plot.pairs, top = text_grob(paste0(gene_name, "\n","p-value = ", p_value, "; log FC = ", logfc), face = "bold", size = 12))
        SaveFigure(plot.pairs,paste0(snp_to_plot,"_",gene_name,"_",i), width = 5, height = 5, res =300)
}

# define a function to plot negative control gRNA-gene pairs
plot_ntc_pair <- function(sceptre_object, response_id, grna_target, p_value) {
	functs_called <- sceptre_object@functs_called
	singleton_integration_strategy <- sceptre_object@grna_integration_strategy == "singleton"
	response_matrix <- get_response_matrix(sceptre_object)
	control_group_complement <- sceptre_object@control_group_complement
	cells_in_use <- sceptre_object@cells_in_use
	# get negative control gRNA index
	grna_group_idxs <- sceptre_object@grna_assignments$indiv_nt_grna_idxs
	grna_group_names <- names(grna_group_idxs)
	y <- sceptre:::load_row(mat = response_matrix, id = response_id)[cells_in_use]
	response_n_umis <- sceptre_object@covariate_data_frame$response_n_umis[cells_in_use]
	normalized_counts <- log(10000 * y/response_n_umis + 1)
	if (length(grna_target == 3)) {
	# get treatment cells contain three negative gRNAs
	nt_idxs <- unique(c(grna_group_idxs[[grna_target[1]]],grna_group_idxs[[grna_target[2]]],grna_group_idxs[[grna_target[3]]]))
	trt_cells <- normalized_counts[nt_idxs]
	} else {
	trt_cells <- normalized_counts[grna_group_idxs[[grna_target]]]
	}
	# get control cells
        cntrl_cells <- normalized_counts[-nt_idxs]
	to_plot <- data.frame(normalized_count = c(trt_cells, cntrl_cells),treatment = factor(c(rep("Treatment", length(trt_cells)), rep("Control", length(cntrl_cells))), levels = c("Treatment", "Control")))
	p_val <- as.numeric(p_value)
	annotation <- as.character(cut(p_val, breaks = c(2, 1, 10^(-seq(2,10)), 0), labels = c(paste0("p <= 10^{", seq(-10, -2),"}"), "p > 0.01", "")))
	set.seed(4)
	to_plot_downsample <- dplyr::sample_n(dplyr::group_by(dplyr::mutate(to_plot,is_zero = (normalized_count == 0)), is_zero, treatment), size = min(dplyr::n(), 1000))
	p_out <- ggplot2::ggplot(data = to_plot, mapping = ggplot2::aes(x = treatment,y = normalized_count, col = treatment)) + ggplot2::geom_violin(draw_quantiles = 0.5,linewidth = 0.6) + ggplot2::geom_jitter(data = to_plot_downsample, alpha = 0.1, size = 0.5) + sceptre:::get_my_theme() + ggplot2::theme(legend.position = "none") + ggplot2::scale_color_manual(values = c("dodgerblue4","firebrick4")) + ggplot2::xlab("Treatment status") + ggplot2::ylab("Normalized expression") + ggplot2::annotate("text",x = 1.5, y = max(to_plot$normalized_count) + 0.5, label = annotation, parse = TRUE) + ggplot2::scale_y_continuous(expand = c(0, 0.1), limits = c(-0.01, max(to_plot$normalized_count) +0.7)) + ggplot2::ggtitle(paste0("Response: ", response_id,"\ngRNA", if (singleton_integration_strategy) "" else " target", ": ", paste0(grna_target[1],"&",grna_target[2],"&",grna_target[3]))) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
	return(p_out)
}


# plot negative control pairs
df_to_plot <- negative_control_pairs[1:15]
for (i in 1:nrow(df_to_plot)) {
	response_id_to_plot <- as.character(df_to_plot[i,"response_id"])
	grna_target_to_plot <- strsplit(as.character(df_to_plot[i,"grna_target"]), "&")[[1]]
	gene_name <- as.character(df_to_plot[i,"gene_symbol"])
	p_value <- sprintf(df_to_plot[i,"p_value"],fmt = "%0.02e")
	logfc = sprintf(df_to_plot[i,"log_2_fold_change"],fmt = "%.2f")
	plot_pairs <- plot_ntc_pair(sceptre_object, response_id_to_plot, grna_target_to_plot, p_value) 
	plot.pairs <- annotate_figure(plot.pairs, top = text_grob(paste0(gene_name, "\n","p-value = ", p_value, "; log FC = ", logfc), face = "bold", size = 12))
        SaveFigure(plot.pairs,paste0("negative_pairs.p",i), width = 5, height = 5, res =300)
}

