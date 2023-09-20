# author: zhaocaiyao


library(DOSE)
library(org.At.tair.db)   ##
library(topGO)
library(clusterProfiler)
#library(pathview)
library(gridExtra)
library(patchwork) 
library(ggplot2)


keytypes(org.At.tair.db) 

args <-commandArgs(T)
#<marker_gene_list><odir><prefix><>
mdate <- gsub("-","",Sys.Date())
prex <- args[3]
prex <- paste(mdate,prex,sep="_")

go_enrich <- function(gene,prefix){
  
  genelits = bitr(gene, 
                  fromType="TAIR", 
                  toType="ENTREZID", 
                  OrgDb="org.At.tair.db")
  # fromType -> toType
  id <- as.vector(genelits[,2])
  #GO analyze
  ego.all <- enrichGO(
    gene          = id,
    keyType = "ENTREZID",
    OrgDb         = org.At.tair.db,
    ont           = "BP",    # optional : "CC" "BP" "MF" "ALL
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE)
  
  # #barplot

p2<-dotplot(ego.all,title=paste("GO Enrichment for ",prefix,sep=""), orderBy = "x")
print(p2)

   return(ego.all@result)
}

setwd(args[2])

genes <- as.vector(read.csv(args[1],header = TRUE,sep="\t"))

mdate<-gsub("-","",Sys.Date())

pdf(paste(mdate,args[3],"Go_clusters.pdf",sep="_"))

for (i in c("W123","W4","W56")){

  gene <- as.vector(genes[which(genes$Group==i),1])
  
  ego_result <- go_enrich(gene, i)
  
  ego_result$Cluster <- rep(i, nrow(ego_result))
  
  write.table(ego_result,file=paste(prex,"go_enrich_stages.csv",sep="_"),append = T,row.names = F, col.names = F,sep = "\t", quote = F)
  
}

dev.off()





