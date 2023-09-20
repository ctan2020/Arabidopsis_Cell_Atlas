
# markerGene_dotplot.R

###############################################################

library(stringr)
library(ggplot2)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(getopt)
library(cowplot)
library(patchwork)
library(readxl)


mydotplot <- function (sample_id,rds_file,work_dir){
    
    setwd(work_dir)

    marker_gene <-  read_xlsx("markerGenes.xlsx",sheet=sample_id) %>% as.data.frame() %>% arrange(Cell_Type,Cluster)

    pre <- Sys.Date() %>% str_replace_all("-","") %>% str_c("_",sample_id,sep="")

    rds <- readRDS(rds_file)
    #csv <- read.csv(marker_gene,header=T)
    csv <- marker_gene 

    p <- DotPlot(rds,features=unique(csv$Gene_ID)) + theme(axis.text.x = element_text(angle = 45,
                                                                               vjust = 0.5, hjust=0.5))
    #match_index <- match(csv$Gene_ID, levels(p$data$features.plot))

    match_index <- match(unique(p$data$features.plot),csv$Gene_ID)

    p$data$features.plot <- factor(p$data$features.plot,levels=unique(p$data$features.plot),labels=csv$Gene_Name[match_index])

    # sorting clusters according to the cell type order
    p$data$id <- factor(p$data$id, levels = unique(csv$Cluster))

    pdfn <- str_c(pre,"_markerGene-dotplot.pdf",sep="")
    #ggsave(pdfn,p,width = 10, height = 6)
    ggsave(pdfn,p,width=16,height =8)
}

rdss <- read_xlsx("markerGenes.xlsx",sheet="rds") %>% as.data.frame()
#i = 1
wd <- "../workdirectory/"
# mydotplot(rdss[i,1],rdss[i,2],wd)


sps <- c(1,2,3,4,5)
#sps <- c(4,14)
for (i in sps){
   mydotplot(sample_id=rdss[i,1],rds_file = rdss[i,2],work_dir = wd)
}







