library(DropletUtils)
library(SoupX)
library(Seurat)

options(future.globals.maxSize = 100000 * 1024^3)
setwd("./")

toc <- Read10X("../filtered_gene_expression",gene.column=1)
tod <- Read10X("../raw_gene_expression",gene.column=1)

#tod <- tod[rownames(toc),]

all <- toc
all <- CreateSeuratObject(all)
#all <- subset(all, subset = nFeature_RNA > 800) 
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)

tod <- tod[rownames(all),]


all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:30)
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = 1:30)

matx <- all@meta.data
toc <- all@assays$RNA@counts
tod <- tod[rownames(all),]
sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc)


out = adjustCounts(sc)
saveRDS(sc,"sc.rds")

#write.csv(out, "matrix.csv")
DropletUtils:::write10xCounts("./soupX_matrix", out,version="3")





