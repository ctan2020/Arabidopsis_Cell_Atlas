library(dplyr)
library(Seurat)

data.list <- list()

##################################################################################################################
s1 <- readRDS('/path/Stage_1_rosette/Stage_1_rosette_harmony_D0_try2_soupx_gene500_10000_lambda1_r0.5/data.cc.RDS')
s1@meta.data$annotation <- as.vector(s1@meta.data$seurat_clusters)
s1@meta.data$annotation[which(s1@meta.data$annotation %in% c(1,4))] = 'Adaxial_PC'
s1@meta.data$annotation[which(s1@meta.data$annotation == 7)] = 'CC'
s1@meta.data$annotation[which(s1@meta.data$annotation == 8)] = 'DC'
s1@meta.data$annotation[which(s1@meta.data$annotation == 2)] = 'Epidermis'
s1@meta.data$annotation[which(s1@meta.data$annotation == 9)] = 'GC'
s1@meta.data$annotation[which(s1@meta.data$annotation %in% c(0,3,5,6,13,14))] = 'Mesophyll'
s1@meta.data$annotation[which(s1@meta.data$annotation == 11)] = 'PC'
s1@meta.data$annotation[which(s1@meta.data$annotation == 10)] = 'PP'
s1@meta.data$annotation[which(s1@meta.data$annotation == 12)] = 'xylem_parencyma'
s1@meta.data$stage <- 'S1'
s1@meta.data$split.clusters <- paste('Stage_1_rosette', s1@meta.data$seurat_clusters)
data.list[[1]] <- s1

##################################################################################################################
s2 <- readRDS('/path/Stage_2_rosette/Stage_2_rosette_harmony_D0_try2_soupx_gene500_10000_lambda1_r0.5/data.cc.RDS')

s2@meta.data$annotation <- as.vector(s2@meta.data$seurat_clusters)
s2@meta.data$annotation[which(s2@meta.data$annotation == 10)] = 'Abaxial_PC'
s2@meta.data$annotation[which(s2@meta.data$annotation == 7)] = 'CC'
s2@meta.data$annotation[which(s2@meta.data$annotation == 8)] = 'DC'
s2@meta.data$annotation[which(s2@meta.data$annotation %in% c(2,3,6))] = 'Epidermis'
s2@meta.data$annotation[which(s2@meta.data$annotation == 9)] = 'GC'
s2@meta.data$annotation[which(s2@meta.data$annotation %in% c(0,1,4))] = 'Mesophyll'
s2@meta.data$annotation[which(s2@meta.data$annotation == 5)] = 'PP'

s2@meta.data$stage <- 'S2'
s2@meta.data$split.clusters <- paste('Stage_2_rosette', s2@meta.data$seurat_clusters)
data.list[[2]] <- s2

##################################################################################################################
s3 <- readRDS('/path/Stage_3_rosette/harmony_D0_try2_soupx_gene500_10000_lambda1_r0.5/data.cc.RDS')
s3@meta.data$annotation <- as.vector(s3@meta.data$seurat_clusters)

s3@meta.data$annotation[which(s3@meta.data$annotation %in% c(4,5))] = 'BS'
s3@meta.data$annotation[which(s3@meta.data$annotation == 7)] = 'CC'
s3@meta.data$annotation[which(s3@meta.data$annotation %in% c(9,13))] = 'DC'
s3@meta.data$annotation[which(s3@meta.data$annotation == 1)] = 'Epidermis'
s3@meta.data$annotation[which(s3@meta.data$annotation == 11)] = 'GC'
s3@meta.data$annotation[which(s3@meta.data$annotation %in% c(0,2,3))] = 'Mesophyll'
s3@meta.data$annotation[which(s3@meta.data$annotation == 12)] = 'PC'
s3@meta.data$annotation[which(s3@meta.data$annotation %in% c(6,10))] = 'PP'
s3@meta.data$annotation[which(s3@meta.data$annotation == 8)] = 'xylem'

s3@meta.data$stage <- 'S3'
s3@meta.data$split.clusters <- paste('Stage_3_rosette', s3@meta.data$seurat_clusters)
data.list[[3]] <- s3

##################################################################################################################
s4 <- readRDS('/path/Stage_4_rosette/Stage_4_rosette_2/harmony_D0_try2_soupx_gene500_10000_lambda1_r0.5/data.cc.RDS')
s4@meta.data$annotation <- as.vector(s4@meta.data$seurat_clusters)
s4@meta.data$annotation[which(s4@meta.data$annotation == 4)] = 'Abaxial_PC'
s4@meta.data$annotation[which(s4@meta.data$annotation == 9)] = 'CC'
s4@meta.data$annotation[which(s4@meta.data$annotation %in% c(1,3))] = 'Epidermis'
s4@meta.data$annotation[which(s4@meta.data$annotation == 12)] = 'GC'
s4@meta.data$annotation[which(s4@meta.data$annotation %in% c(0,2,7,8))] = 'Mesophyll'
s4@meta.data$annotation[which(s4@meta.data$annotation == 10)] = 'PP'
s4@meta.data$annotation[which(s4@meta.data$annotation == 6)] = 'Vascular'
s4@meta.data$annotation[which(s4@meta.data$annotation == 10)] = '5'

s4@meta.data$stage <- 'S4'
s4@meta.data$split.clusters <- paste('Stage_4_rosette', s4@meta.data$seurat_clusters)
data.list[[4]] <- s4
##################################################################################################################
S4_leaf_base	<- readRDS('/path/S4_leaf_segment/S4_leaf_base/harmony_S4_leaf_try1_2samples_gene500_10000_lambda1_r0.5/data.cc.RDS')
S4_leaf_base@meta.data$annotation <- as.vector(S4_leaf_base@meta.data$seurat_clusters)
S4_leaf_base@meta.data$annotation[which(S4_leaf_base@meta.data$annotation == 3 )] = 'Epidermis'
S4_leaf_base@meta.data$annotation[which(S4_leaf_base@meta.data$annotation %in% c(0,1,4,7,8,9))] = 'Mesophyll'
S4_leaf_base@meta.data$annotation[which(S4_leaf_base@meta.data$annotation == 2)] = 'Vascular'
S4_leaf_base@meta.data$stage <- 'S4_leaf_base'
S4_leaf_base@meta.data$split.clusters <- paste('S4_leaf_base', S4_leaf_base@meta.data$seurat_clusters)
data.list[[7]] <- S4_leaf_base
##################################################################################################################
S4_leaf_middle	<- readRDS('/path/S4_leaf_segment/S4_leaf_middle/harmony_S4_leaf_try1_2samples_gene500_10000_lambda1_r0.5/data.cc.RDS')
S4_leaf_middle@meta.data$annotation <- as.vector(S4_leaf_middle@meta.data$seurat_clusters)
S4_leaf_middle@meta.data$annotation[which(S4_leaf_middle@meta.data$annotation == 0 )] = 'Epidermis'
S4_leaf_middle@meta.data$annotation[which(S4_leaf_middle@meta.data$annotation %in% c(1,2,3,5,7))] = 'Mesophyll'
S4_leaf_middle@meta.data$annotation[which(S4_leaf_middle@meta.data$annotation %in% c(4,9))] = 'Vascular'
S4_leaf_middle@meta.data$stage <- 'S4_leaf_middle'
S4_leaf_middle@meta.data$split.clusters <- paste(' S4_leaf_middle', S4_leaf_middle@meta.data$seurat_clusters)
data.list[[8]] <- S4_leaf_middle
##################################################################################################################
S4_leaf_tip	<- readRDS('/path/S4_leaf_segment/S4_leaf_tip/harmony_S4_leaf_try1_2samples_gene500_10000_lambda1_r0.5/data.cc.RDS')
S4_leaf_tip@meta.data$annotation <- as.vector(S4_leaf_tip@meta.data$seurat_clusters)
S4_leaf_tip@meta.data$annotation[which(S4_leaf_tip@meta.data$annotation %in% c(0,2,6))] = 'Epidermis'
S4_leaf_tip@meta.data$annotation[which(S4_leaf_tip@meta.data$annotation %in% c(1,3,4,9,10))] = 'Mesophyll'
S4_leaf_tip@meta.data$annotation[which(S4_leaf_tip@meta.data$annotation == 11)] = 'Vascular'
S4_leaf_tip@meta.data$stage <- 'S4_leaf_tip'
S4_leaf_tip@meta.data$split.clusters <- paste('S4_leaf_tip', S4_leaf_tip@meta.data$seurat_clusters)
data.list[[9]] <- S4_leaf_tip
##################################################################################################################
s5 <- readRDS('/path/Stage_5_rosette/Stage_5_rosette_harmony_D0_try2_soupx_gene500_10000_lambda1_r0.5/data.cc.RDS')
s5@meta.data$annotation <- as.vector(s5@meta.data$seurat_clusters)

s5@meta.data$annotation[which(s5@meta.data$annotation %in% c(1))] = 'BS'
s5@meta.data$annotation[which(s5@meta.data$annotation %in% c(7))] = 'BS-PC'
s5@meta.data$annotation[which(s5@meta.data$annotation %in% c(6))] = 'BS-XP'
s5@meta.data$annotation[which(s5@meta.data$annotation %in% c(8))] = 'CC'
s5@meta.data$annotation[which(s5@meta.data$annotation %in% c(2))] = 'Epidermis'
s5@meta.data$annotation[which(s5@meta.data$annotation %in% c(0,3,4))] = 'Mesophyll'
s5@meta.data$annotation[which(s5@meta.data$annotation %in% c(5))] = 'PP'
s5@meta.data$annotation[which(s5@meta.data$annotation %in% c(9))] = 'Vascular'

s5@meta.data$stage <- 'S5'
s5@meta.data$split.clusters <- paste('Stage_5_rosette', s5@meta.data$seurat_clusters)
data.list[[5]] <- s5

##################################################################################################################
s6 <- readRDS('/path/Stage_6_rosette/harmony_D0_try2_soupx_gene500_10000_lambda1_r0.5/data.cc.RDS')
s6@meta.data$annotation <- as.vector(s6@meta.data$seurat_clusters)
s6@meta.data$annotation[which(s6@meta.data$annotation %in% c(5,7,8))] = 'BS'
s6@meta.data$annotation[which(s6@meta.data$annotation %in% c(4))] = 'BS-PC'
s6@meta.data$annotation[which(s6@meta.data$annotation %in% c(9,11))] = 'CC'
s6@meta.data$annotation[which(s6@meta.data$annotation %in% c(1))] = 'Epidermis'
s6@meta.data$annotation[which(s6@meta.data$annotation %in% c(0,2))] = 'Mesophyll'
s6@meta.data$annotation[which(s6@meta.data$annotation %in% c(12))] = 'PC'
s6@meta.data$annotation[which(s6@meta.data$annotation %in% c(10))] = 'PP'
s6@meta.data$annotation[which(s6@meta.data$annotation %in% c(6))] = 'Vascular'
s6@meta.data$annotation[which(s6@meta.data$annotation %in% c(9))] = 'CC'


s6@meta.data$stage <- 'S6'
s6@meta.data$split.clusters <- paste('Stage_6_rosette', s6@meta.data$seurat_clusters)
data.list[[6]] <- s6

adata <- merge(data.list[[1]], y=data.list[-1])
adata <- NormalizeData(object = adata,normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(adata,'rosette_combine_0710.rds')


adata$a_s<-paste0(adata$annotation,'_',adata$stage)
exp <- AverageExpression(adata,group.by = "a_s")[['RNA']]
write.table(exp,'rosette_exp_0710.txt',sep="\t",quote=FALSE,row.names=TRUE)








