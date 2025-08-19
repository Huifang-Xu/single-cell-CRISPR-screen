library(ggplot2)
library(data.table)
library(dplyr)

outputDir <- "/scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/batch2" 

guide_assignment <- fread("/scratch/hx37930/project/MTAG/06.crispr/results/scRNA_crispri_merged/combined_crispr/yk3/guide_RNAs_filtered/guide_assignment.csv")

n_grna_cell <- data.frame(table(guide_assignment$bc_wells))

# calculate cumulative percentage
n_grna_cell.pct <- n_grna_cell %>%
  group_by(Freq) %>%
  summarise(proportion = n()) %>%
  mutate(Perc = cumsum(100*proportion/sum(proportion))) %>% select(-proportion)

png(paste0(outputDir,"/dist_gRNA.v1.png"),type=c("cairo", "cairo-png", "Xlib", "quartz"), width = 6, height = 4, units = "in", res = 400)
ggplot(n_grna_cell) + geom_histogram(aes(x = Freq),binwidth=1,fill = "darkseagreen")+theme_bw()+xlab("# gRNAs per cell")+ ylab("# cells")
dev.off()

png(paste0(outputDir,"/dist_gRNA.v2.png"),type=c("cairo", "cairo-png", "Xlib", "quartz"), width = 6, height = 4, units = "in", res = 400)
ggplot(n_grna_cell) + geom_histogram(aes(x = Freq),binwidth=1,fill = "darkseagreen")+theme_bw()+xlab("# gRNAs per cell")+ ylab("# cells")+xlim(0,32)
dev.off()

png(paste0(outputDir,"/dist_gRNA_perc.v1.png"),type=c("cairo", "cairo-png", "Xlib", "quartz"), width = 6, height = 4, units = "in", res = 400)
ggplot(n_grna_cell.pct)+ geom_point(aes(x =Freq,y=Perc))+theme_bw()+xlab("# gRNAs per cell")+ ylab("% cells")
dev.off()

png(paste0(outputDir,"/dist_gRNA_perc.v2.png"),type=c("cairo", "cairo-png", "Xlib", "quartz"), width = 6, height = 4, units = "in", res = 400)
ggplot(n_grna_cell.pct)+ geom_point(aes(x =Freq,y=Perc))+theme_bw()+xlab("# gRNAs per cell")+ ylab("% cells")+xlim(0,50)
dev.off()
