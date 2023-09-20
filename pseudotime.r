library(Seurat)
library(monocle)
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(scales)
library(RColorBrewer)

args <- commandArgs(T)
# [1-sample.rds 2-markergene.list]
ATNL2729 <- readRDS(args[1])
data.all <- as.matrix(ATNL2729@assays[["RNA"]]@data)
meta.epi <- as.matrix(ATNL2729@meta.data[["seurat_clusters"]])
rownames(meta.epi) <- colnames(data.all)
meta.epi <- as.matrix(meta.epi)
data_matrix <- as.matrix(data.all[,rownames(meta.epi)])
feature_ann <- data.frame(gene_id=rownames(data_matrix),gene_short_name=rownames(data_matrix))
rownames(feature_ann) <- rownames(data_matrix)
data_fd <- new("AnnotatedDataFrame", data = feature_ann)

sample_ann <- meta.epi
sample_ann <- as.data.frame(sample_ann)
data_pd <- new("AnnotatedDataFrame", data =sample_ann)

data.cds <- newCellDataSet(as(data_matrix,'sparseMatrix'),phenoData =data_pd,featureData =data_fd,expressionFamily=negbinomial.size())
data.cds <- estimateSizeFactors(data.cds)
data.cds <- estimateDispersions(data.cds)

expressed_genes <- row.names(subset(fData(data.cds)))
data.cds@phenoData@data$V1 <- factor(data.cds@phenoData@data$V1 )

ordering_genes <- read_tsv(argv[2])
ordering_genes <- ordering_genes$Gene
data.cds <- setOrderingFilter(data.cds, ordering_genes)
#plot_ordering_genes(data.cds)

data.cds <- reduceDimension(data.cds,max_components = 2, method = 'DDRTree',max_components = 2, norm_method = 'log')
data.cds <- orderCells(data.cds,reverse=F)

#color plot
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

pdata <- read.csv('pData2.1.csv')
data.cds@phenoData@data$Time <- pdata$Time
data.cds@phenoData@data$Time <- as.factor(data.cds@phenoData@data$Time)
data.cds@phenoData@data$RawCluster <- pdata$Rcluster
data.cds@phenoData@data$RawCluster <- as.factor(data.cds@phenoData@data$RawCluster)

pdf("Pseudotime_pollen.times.1.pdf",width = 8, height = 5)
plot(plot_cell_trajectory(data.cds, show_cell_names = F, color_by = "Pseudotime")+scale_color_viridis_c())
plot(plot_cell_trajectory(data.cds, show_cell_names = F, color_by = "V1")+scale_color_manual(values = getPalette(length(unique(data.cds@phenoData@data[,"V1"])))))
plot(plot_cell_trajectory(data.cds, show_cell_names = F, color_by = "State")+scale_color_manual(values = getPalette(length(unique(data.cds@phenoData@data[,"State"])))))
plot(plot_cell_trajectory(data.cds, show_cell_names = F, color_by = "Time") +scale_color_manual(values = getPalette(length(unique(data.cds@phenoData@data[,"Time"])))))
plot(plot_cell_trajectory(data.cds, show_cell_names = F, color_by = "RawCluster") +scale_color_manual(values = getPalette(length(unique(data.cds@phenoData@data[,"RawCluster"])))))
dev.off()


##plot markergenes
diff_test_res <- differentialGeneTest(data.cds[expressed_genes,],fullModelFormulaStr ="~sm.ns(Pseudotime)")
diff_test_res <- diff_test_res[order(diff_test_res$qval),]
write.csv(diff_test_res,file = "pollen_pseudotime_markergenes.csv")


# Clustering top 100 genes by pseudotemporal expression pattern
pdf("Pseudotime_pollen.times.3.pdf")
if(nrow(diff_test_res)>100){
sig_gene_names <- row.names(diff_test_res[1:100,]) # Select top 100 gene used to cluster
}else{
sig_gene_names <- row.names(diff_test_res)
}
plot_pseudotime_heatmap(data.cds[sig_gene_names,],num_clusters = 3, cores = 1,show_rownames = T)
dev.off()

pdf("Pseudotime_pollen.times.4.pdf",width = 12 , height = 40)
if(nrow(diff_test_res)>200){
sig_gene_names <- row.names(diff_test_res[1:200,])
}else{
sig_gene_names <- row.names(diff_test_res)
}
plot_pseudotime_heatmap(data.cds[sig_gene_names,],num_clusters = 3, cores = 1,show_rownames = T)
dev.off()

pdf("Pseudotime_pollen.times.5.pdf",width = 12 , height = 40)
if(nrow(diff_test_res)>500){
sig_gene_names <- row.names(diff_test_res[1:500,])
}else{
sig_gene_names <- row.names(diff_test_res)
}
plot_pseudotime_heatmap(data.cds[sig_gene_names,],num_clusters = 3, cores = 1,show_rownames = T)
dev.off()

pdf("Pseudotime_pollen.times.2.2.pdf")
cds_subset <- data.cds[rownames(diff_test_res)[1:20],]
plot_genes_in_pseudotime(cds_subset, ncol = 4, color_by = 'Time')+scale_color_manual(values = getPalette(length(unique(data.cds@phenoData@data[,'Time']))))
dev.off()











