cat>"sc_Seurat_arth.R"<<EOF

library(getopt)
arg <- matrix(c("input", "i","1","character","input file1",
                "outdir","o","1","character","outdir",
                "sample","s","1","character","sample,default=Maize",
                "tissue","t","1","character","tissue,default=Embro",
                "cellcycle","c","1","character","either remove the effect of cell cycle genes or not",
		"regressall","a","1","character","regress all the effects or half",
		"cc.genes","e","1","character","input cell cycle genes list",
		"dims","d","1","integer","dims option for FindNeighbors,default=15",
                "resolution","r","1","numeric","resolution option for FindClusters,[0.4-1.2],default=1",
                "help","h","0","logical", "Usage: Rscript runiDrop.R -i <input> -o <outdir> [-s SAMPLE -t TISSUE]",
                "minCG","m","1","integer","minimum gene number for each cell, i.e. nFeature_RNA, default=500",
                "maxCG","n","1","integer"," max gene number for each cell, i.e. nFeature_RNA, default=5000",
                "thedata","f","0","logical", "Save the RDS file",
                "lambda","l","1","numeric"," RUNharmony lambda, default=1",
                "featurePlot","p","1","character","T or F",
                "marker.colname","k","1","character","markergene colnames",
		"features","g","1","character","plot genes expression"
		),byrow=T,ncol=5)
opt = getopt(arg)
if (is.null(opt$sample)){
        opt$sample <- "Maize"
}
if (is.null(opt$tissue)){
        opt$tissue <- "stem"
}
if (is.null(opt$outdir)){
        opt$outdir <- "output"
}
if (is.null(opt$cellcycle)){
        opt$cellcycle <- F
}
if (is.null(opt$regressall)){
        opt$regressall <- T
}

if (is.null(opt$dims)){
        opt$dims <- 15
}
if (is.null(opt$resolution)){
        opt$resolution <- 1
}
if (is.null(opt$minCG)){
    opt$minCG <- 500
}
if (is.null(opt$maxCG)){
    opt$maxCG <- 5000
}


if (is.null(opt$lambda)){
    opt$lambda <- 1
}
if (is.null(opt$featurePlot)){
 opt$featurePlot <- F
}

library(harmony)
library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(Matrix)
library(DoubletFinder)

file=scan(opt$input,"")
#print(file)

##make function of doublefinkder 
Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(0.05*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  #data<-subset(data,subset=doublet_info=="Singlet")
  data
}


#dir.create(opt$outdir);
setwd(opt$outdir);
thedatalist<-list()
n=0
w=c()
h=c()
id=c()
vf<-c()
plot.list<-list()
thedata <- c()
for(i in file)
{
n=n+1
matrix_dir=i
#setwd(matrix_dir)
print(matrix_dir)
mat <- Read10X(matrix_dir,gene.column=1)
thedata <- CreateSeuratObject(counts = mat, assay = "RNA",min.cells = 3, project = opt$sample, min.features = opt$minCG)
print(thedata)
thedata[["percent.atm"]] <- PercentageFeatureSet(thedata, pattern = "^ATM")
thedata[["percent.atc"]] <- PercentageFeatureSet(thedata, pattern = "^ATC")

    thedata <- subset(thedata, subset = nFeature_RNA > opt$minCG & nFeature_RNA < opt$maxCG & percent.atm < 5 & percent.atc <5)
print(thedata)
    name1 <- strsplit(i,"/",fixed=T)[[1]][length(strsplit(i,"/",fixed=T)[[1]])-1]
	#name1 <- B73_L3_web_0
#    name2 <- gsub("_\\d.*web_.*","",name1)
	thedata@meta.data$batch <- name1
#	thedata@meta.data$group <- name2
    #thedata[[n]]@meta.data$section<-opt$section
    #thedata[[n]]@meta.data$DAP_list<-DAP_list[[n]]
     thedata <- NormalizeData(thedata)
     thedata <- FindVariableFeatures(thedata, selection.method = "vst", nfeatures = 2000)
     thedata <- ScaleData(thedata)
     thedata <- RunPCA(thedata)
     thedata <- RunUMAP(thedata, dims = 1:30)

	thedata<-Find_doublet(thedata)
	thedata<-subset(thedata,subset=doublet_info=="Singlet")


    thedatalist[[n]] <- thedata
    id[n]<- name1
    plot1 <- VlnPlot(thedata, features = c("nFeature_RNA", "nCount_RNA"))
    plot.list[[n]] <- plot1
}        
#saveRDS(thedata,"data_list.RDS")

pdf("01.genes_UMIcounts.plot.pdf")
plot.list
dev.off()


data.harmony<-merge(thedatalist[[1]], y=thedatalist[-1],add.cell.ids = id)
data.harmony <- NormalizeData(data.harmony, normalization.method = "LogNormalize", scale.factor = 10000)
data.harmony <- FindVariableFeatures(data.harmony,selection.method = "vst", nfeatures = 3000,verbose = TRUE,)
data.harmony <- ScaleData(data.harmony, features = rownames(data.harmony),verbose = TRUE)
data.harmony <- RunPCA(data.harmony, npcs = 30, verbose = FALSE)

#SCT or RNA
set.seed(100)
system.time({data.harmony <- RunHarmony(data.harmony, 'batch', plot_convergenc=TRUE,lambda = opt$lambda)})

pdf("02_ElbowPlot.pdf", width = 10)
ElbowPlot(object = data.harmony)
dev.off()

set.seed(100)
data.harmony <- FindNeighbors(data.harmony, reduction = "harmony",dims = 1:opt$dims)
data.harmony <- FindClusters(data.harmony, resolution = opt$resolution)
data.harmony <- RunUMAP(data.harmony, reduction = "harmony",dims = 1:opt$dims)
data.harmony <- RunTSNE(data.harmony, reduction = "harmony",dims = 1:opt$dims)


pdf("03_umap.tsne.pdf", width = 24, height=8)
plot1 <- DimPlot(object = data.harmony, reduction = "umap", label = T,pt.size=0.5, raster = FALSE, group.by = c("ident", "batch"))
plot1.1 <- FeaturePlot(data.harmony, features = "nFeature_RNA", reduction = "umap", label = T, combine = F, order = T)
plot1+plot1.1
plot2 <- DimPlot(object = data.harmony, reduction = "tsne", label = T,pt.size=0.5, raster = FALSE, group.by = c("ident", "batch"))
plot2.1 <- FeaturePlot(data.harmony, features = "nFeature_RNA", reduction = "tsne", label = T, combine = F, order = T)
plot2+plot2.1
dev.off()


pdf("04_umap.tsnesplit.pdf", width = 24, height=6)
plot3 <- DimPlot(object = data.harmony, reduction = "umap", pt.size=1,group.by = "ident",split.by = "batch")
plot3
plot4 <- DimPlot(object = data.harmony, reduction = "tsne", pt.size=1,group.by = "ident",split.by = "batch")
plot4
dev.off()

saveRDS(data.harmony,"data.RDS")

if(opt$cellcycle){
        #input cell cycle genes
        cc.genes <- read.csv(opt$cc.genes, header = T)
        s.genes <- cc.genes[which(cc.genes[,1]=="S"),"Gene"]
        s.genes
        g2m.genes <- c(cc.genes[which(cc.genes[,1] =="G2"),"Gene"], cc.genes[which(cc.genes[,1] =="M"),"Gene"])
        g2m.genes

        ##Assign Cell-Cycle Scores
        data.harmony.cc <- CellCycleScoring(data.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

        #view cell cycle scores and phase assignments
        head(data.harmony.cc[[]])

        if(opt$regressall){
        ##Regress out cell cycle scores during data scaling
        data.harmony.cc <- ScaleData(data.harmony.cc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(data.harmony.cc))
        } else {
          data.harmony.cc$CC.Difference <- data.harmony.cc$S.Score - data.harmony.cc$G2M.Score
          data.harmony.cc <- ScaleData(data.harmony.cc, vars.to.regress = "CC.Difference", features = rownames(data.harmony.cc))
        }
	data.harmony.cc <- RunPCA(data.harmony.cc, features = VariableFeatures(data.harmony), nfeatures.print = 10)
        # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
       #data.harmony <- RunPCA(data.harmony, features = c(s.genes, g2m.genes))
        #DimPlot(data.harmony)
	
        ###reharmony
        data.harmony.cc <- RunHarmony(data.harmony.cc, 'batch', plot_convergenc=TRUE,lambda = opt$lambda)
        data.harmony.cc <- RunUMAP(data.harmony.cc, reduction = "harmony", dims = 1:opt$dims)
        data.harmony.cc <- FindNeighbors(data.harmony.cc, reduction = "harmony", dims = 1:opt$dims)
        data.harmony.cc <- FindClusters(data.harmony.cc, resolution = opt$resolution)

        plot3 <- DimPlot(data.harmony.cc, reduction = "umap", raster = FALSE, group.by = c("Phase","seurat_clusters","old.ident"), label = TRUE, repel = TRUE)
	plot2 <- DimPlot(object = data.harmony.cc, reduction = "umap", label = T, pt.size=0.5, raster = FALSE, group.by = c("ident",  "batch"))
	plot2.1 <- FeaturePlot(data.harmony.cc, features = "nFeature_RNA", reduction = "umap", label = T, combine = F, order = T)

        pdf("03.cc_umap.harmony.pdf", height = 8, width = 24)
        print(plot3)
	print(plot2+plot2.1)
        dev.off()
		
	pdf("04.cc_umap.harmony.split.pdf", width = 24, height=6)
	plot4 <- DimPlot(object = data.harmony, reduction = "umap", pt.size=1,group.by = "ident",split.by = "batch")
	print(plot4)
	dev.off()
        data.harmony
        saveRDS(data.harmony.cc, "data.cc.RDS")


        pdf("05.cc_cell_number.pdf", width = 12, height = 8)
        dat<-as.data.frame(table(data.harmony.cc@active.ident))
        p<-ggplot(dat,aes(Var1,Freq,fill=Var1))+geom_bar(stat = "identity",position = "dodge") + geom_text(aes(label=Freq,vjust=-0.5) )
        p
        dev.off()

        write.table(dat,file = "05.cc_cell_number.txt",sep="\t",quote=F,row.names = F)

        ###all marker genes
        DefaultAssay(data.harmony.cc) <- "RNA"
        adata.markers <- FindAllMarkers(data.harmony.cc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
        allmakers <- adata.markers %>% group_by(cluster)
        write.csv(allmakers,file = "06.cc.markergenes_list6.csv")

        ###top10 marker genes heatmap
        pdf('06.cc.top10_markergenes.heatmap.integreted.pdf',width = 12 , height = 8 )
        top10 <- adata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
        heatmap <- DoHeatmap(data.harmony.cc, features = top10$gene)
        print(heatmap)
        dev.off()
}
EOF


try=D0_try2_soupx
tissue=D0_silique

for minCG in 500; do
for maxCG in 10000; do
for lambda in 1 ;do


mkdir harmony_${try}_gene${minCG}_${maxCG}_lambda${lambda}_r0.5
cd harmony_${try}_gene${minCG}_${maxCG}_lambda${lambda}_r0.5
cp ../soupx.file ./
cat>"sc_Seurat_harmony.sh"<<EOF
export R_LIBS_USER=/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/wangfang/02.software/miniconda3/envs/R_4.1/lib/R/library
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/wangfang/02.software/miniconda3/envs/R_4.2/bin/Rscript \
/hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhaocaiyao/scRNA_seq/script/new.data.script/sc_Seurat_arth.R \
-i soupx.file -o ./ \
-s Arth -t $tissue \
-c T \
-a T \
-e /hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhaocaiyao/C3_C4/marker_genes/cc.genes.arth.2.maize1.csv \
-d 30 -r 0.5 -m ${minCG} -n ${maxCG} -l ${lambda} \
-p T \
-k marker_gene_ID \
-g /hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhaocaiyao/scRNA_seq/Find.cluster/markergenes/arth.atlas.all.markergenes.csv
EOF

qsub -clear -cwd -q st.q -P P20Z10200N0035 -l vf=200g,num_proc=2 -binding linear:2 sc_Seurat_harmony.sh

cd ../
done
done
done
