library("WGCNA")
library("reshape2")
library("stringr")
library("ggplot2")
options(stringsAsFactors=FALSE)
enableWGCNAThreads(nThreads=10)
data<-read.table("S1-S6.all.exp.new.txt",header=TRUE,row.names=1)
data<-data[which(rowSums(data)>0),]*1.0
type="unsigned"
corType="pearson"
corFnc=ifelse(corType=="pearson", cor, bicor)
maxPOutliers=ifelse(corType=="pearson",1,0.05)
robustY=ifelse(corType=="pearson",T,F)
dim(data)
m.mad <- apply(data,1,mad)
dataExprVar <-data[which(m.mad>max(quantile(m.mad,probs=seq(0,1,0.25))[2],0.01)),]
dataExpr<-as.data.frame(t(dataExprVar))
gsg=goodSamplesGenes(dataExpr, verbose=3)
if (!gsg$allOK){
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", 
                paste(names(dataExpr)[!gsg$goodGenes], collapse=",")));
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", 
                paste(rownames(dataExpr)[!gsg$goodSamples], collapse=",")));
        dataExpr=dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes=ncol(dataExpr)
nSamples=nrow(dataExpr)
dim(dataExpr)
sampleTree=hclust(dist(dataExpr), method="average")
pdf(file="sample_clustering.pdf")
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="")
dev.off()
powers=c(c(1:10), seq(from=12, to=30, by=2))
sft=pickSoftThreshold(dataExpr, powerVector=powers,networkType=type, verbose=5)
pdf(file="Topology_model.pdf")
par(mfrow=c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main=paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,cex=cex1, col="red")
dev.off()
power=sft$powerEstimate
power

if (is.na(power)){
    power=ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
    ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
            ifelse(type == "unsigned", 6, 12))       
        )
    )
}
blockSize(nGenes,rectangularBlocks=TRUE,maxMemoryAllocation=2^30) 
net=blockwiseModules(dataExpr, power=power, maxBlockSize=nGenes,TOMType=type, minModuleSize=50,reassignThreshold=0, mergeCutHeight=0.3,numericLabels=TRUE, pamRespectsDendro=FALSE,saveTOMs=FALSE, corType=corType, maxPOutliers=maxPOutliers, loadTOMs=TRUE,verbose=3)
table(net$colors)
mergedColors=labels2colors(net$colors)
table(mergedColors)

moduleLabels=net$colors
moduleColors=labels2colors(moduleLabels)
pdf(file="module_colors.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels=FALSE, hang=0.03,addGuide=TRUE, guideHang=0.05)
dev.off()
pdf(file="merged_module_colors.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels=FALSE, hang=0.03,addGuide=TRUE, guideHang=0.05)
dev.off()

MEs=net$MEs
MEs_col=MEs
colnames(MEs_col)=paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col=orderMEs(MEs_col)
write.csv(dataExpr,"data.tom")

trait<-read.table(file="S1-S6.SAG-growth.index.info",header=T,row.names=1)
MEs_colpheno=orderMEs(cbind(MEs_col,trait))
pdf(file="Eigengene-trait_adj.heatmap.pdf",width=10,height=10)
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", marDendro=c(4,4,4,4),marHeatmap=c(4,4,4,4), plotDendrograms=T,xLabelsAngle=90)
dev.off()

adjacency=adjacency(dataExpr,power=power)
TOM=TOMsimilarity(adjacency)
dissTOM=1-TOM
geneTree=hclust(as.dist(dissTOM),method="average")
dynamicMods=cutreeDynamic(dendro=geneTree,distM=dissTOM,deepSplit=2,pamRespectsDendro=FALSE,minClusterSize=30)
dynamicColors=labels2colors(dynamicMods)

merge_modules=mergeCloseModules(dataExpr,dynamicColors,cutHeight=0.1,verbose=3)
mergedColors=merge_modules$colors
mergedMEs=merge_modules$newMEs
pdf("module cluster.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
    c("Dynamic Tree Cut", "Merged dynamic"),
    dendroLabels=FALSE, hang=0.03,
    addGuide=TRUE, guideHang=0.05)
dev.off()

dataKME=signedKME(dataExpr,MEs,outputColumnName="kME_MM.")
write.csv(dataKME, "kME_MM.csv")
sag=as.data.frame(trait$sagindex)
genTraitSig2=as.data.frame(cor(dataExpr,sag,use="p"))
GSPvalue2=as.data.frame(corPvalueStudent(as.matrix(genTraitSig2),nSamples))
write.csv(genTraitSig2,"GS.SAGindex.csv")
write.csv(GSPvalue2,"GS.SAGindex.csv")

probes=colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
cyt=exportNetworkToCytoscape(TOM,edgeFile="all.edges.txt",nodeFile="all.nodes.MES.txt",weighted=TRUE,threshold=0,nodeNames=probes,nodeAttr=moduleColors)
