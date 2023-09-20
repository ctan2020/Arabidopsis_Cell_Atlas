rds=/hwfssz1/ST_EARTH/P20Z10200N0035/USER/xiangsunhuan/1.scRNA_seq/01.Arabidopsis/02.analysis/01.clustering/02.V2_data_20230203/6R_1117/6R_1117_harmony_D0_try2_soupx_gene500_10000_lambda1_r0.5/data.cc.RDS
clusters=seurat_clusters
proname=D16_root
cat>average.exp.correlation.R<<EOF
library(Seurat)
library(ggplot2)
rds <- readRDS("$rds")
cell.num <- data.frame(table(rds@meta.data[["$clusters"]]))
write.table(cell.num, paste0("$proname",".cell.num.txt"))
#Idents(rds) <- rds@meta.data[["$clusters"]]
ave.exp <- AverageExpression(rds, group.by = "seurat_clusters", slot = "data")[["RNA"]]
colnames(ave.exp) <- as.character(colnames(ave.exp))
write.table(ave.exp, paste0("$proname", ".average.expression.txt"))
cell.ann <- read.csv("$ann")
cell.ann <- distinct(cell.ann, Cluster, .keep_all = TRUE)
y <- as.character(cell.ann[,"Cluster"])
ave.exp <- ave.exp[,y]
colnames(ave.exp) <- paste(cell.ann[,"Cell_Type"], cell.ann[,"Cluster"], sep = "_")
cor.exp <- as.data.frame(cor(ave.exp))
cor.exp[["x"]] <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, -x)
p1 <- ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile() + scale_fill_gradient(high = "red", low = "white")
p1 <- p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1$"data"$"x" <- forcats::fct_inorder(p1$"data"$"x")
p1$"data"$"y" <- forcats::fct_inorder(p1$"data"$"y")
pdf(paste0("$proname", ".seurat.cluster.correlation.heatmap.pdf"), width = 9.55, height = 8)
p1
dev.off()

EOF

cat>ave.exp.sh<<EOF
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript average.exp.correlation.R
EOF

#sh ave.exp.sh
qsub -clear -cwd -q st.q -P P20Z10200N0035 -l vf=40g,num_proc=1 -binding linear:1 ave.exp.sh

