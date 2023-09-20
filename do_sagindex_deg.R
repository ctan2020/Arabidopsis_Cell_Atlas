

library(Seurat)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
#library(hrbrthemes)
library(viridis)
library(cowplot)
library(openxlsx)

sag_degs6 <- function(args){
	
	irds <- args[1]
	sp <- args[2]

	pre <- Sys.Date() %>% str_replace_all("-","") %>% str_c("_sagvalue-comparison_",sp)

	# step 1: extract top 100 and bottom 100 cells with sag-index value
	
	icell <- "/hwfssz1/ST_EARTH/P20Z10200N0035/USER/tancong/scRNA/1.leaf_development/6.SAG_index/20211011_SAG-index.csv"
	icell <- read_csv(icell)
	# Row.names,Cell_ID,Stage,sagindex,Cell

#	icell1 <- icell %>% filter(Cell==sp) %>% slice_max(order_by=sagindex,n=100) %>% mutate(sag_group="High")
#	icell2 <- icell %>% filter(Cell==sp) %>% slice_min(order_by=sagindex,n=100) %>% mutate(sag_group="Low")
	icell_w3 <- icell %>% filter(Cell==sp & sagindex <5 &  Stage=="W3")
	icell_w6 <- icell %>% filter(Cell==sp & sagindex <5 &  Stage=="W6")
	icell <- icell_w6
	median100_max <- icell$sagindex[icell$sagindex > icell$sagindex %>% median()] %>% sort() %>% head(50) %>% max()
	median100_min <- icell$sagindex[icell$sagindex < icell$sagindex %>% median()] %>% sort() %>% tail(50) %>% min()

	icell1 <- icell %>% filter(Cell==sp) %>% slice_max(order_by=sagindex,n=100) %>% mutate(sag_group="High")
	icell1 %>% summary()
	
	icell2 <- icell %>% filter(sagindex >= median100_min & sagindex <= median100_max ) %>% mutate(sag_group="Low")
	icell2 %>% summary()
		
	icell <- rbind(icell1,icell2) %>% as.data.frame()
	tsvn1 <- str_c(pre,"_cell-sag_W6.csv")
        write_tsv(icell,tsvn1)

	# plot
	p <- ggplot(data=icell,aes(sag_group,sagindex,fill=sag_group)) + geom_boxplot() +  theme_cowplot()+  theme(legend.position="none",plot.title = element_text(size=11)) + ggtitle("boxplot of sag-index value for W6  senescence cells")
	ggsave(str_c(pre,"-W6-","boxplot.png"),p, dpi=600)
		
	
	# step 2: prepare rds
	irds <- readRDS(irds)
	rds <- subset(irds,subset=(Cell_ID %in% unique(icell$Cell_ID )))
	
	# step 3: DEGs detection
	tmp <- rds@meta.data
  
  tmp <- tmp %>% left_join(icell,by=c("Cell_ID"="Cell_ID")) %>% select(sag_group)
	
	rds@meta.data$sag_group <- tmp$sag_group
	
	Idents(rds) <- rds@meta.data$sag_group
	
	markers <- FindAllMarkers(rds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	
	markers <- markers %>% group_by(cluster) %>% as.data.frame()
	
	tsvn <- str_c(pre,"_markers_W6.csv")
	write_tsv(markers,tsvn)

}


sag_degs3 <- function(args){
	
	irds <- args[1]
	sp <- args[2]
	oxlsx <- list()
	pre <- Sys.Date() %>% str_replace_all("-","") %>% str_c("_sagvalue-comparison_",sp)

	# step 1: extract top 100 and bottom 100 cells with sag-index value
	
	icell <- "/hwfssz1/ST_EARTH/P20Z10200N0035/USER/tancong/scRNA/1.leaf_development/6.SAG_index/20211011_SAG-index.csv"
	icell <- read_csv(icell)
	# Row.names,Cell_ID,Stage,sagindex,Cell

#	icell1 <- icell %>% filter(Cell==sp) %>% slice_max(order_by=sagindex,n=100) %>% mutate(sag_group="High")
#	icell2 <- icell %>% filter(Cell==sp) %>% slice_min(order_by=sagindex,n=100) %>% mutate(sag_group="Low")
	
	icell_w3 <- icell %>% filter(Cell==sp & sagindex <5 &  Stage=="W3")
	icell_w6 <- icell %>% filter(Cell==sp & sagindex <5 &  Stage=="W6")
	icell <- icell_w3
	median100_max <- icell$sagindex[icell$sagindex > icell$sagindex %>% median()] %>% sort() %>% head(50) %>% max()
        median100_min <- icell$sagindex[icell$sagindex < icell$sagindex %>% median()] %>% sort() %>% tail(50) %>% min()

        #icell1 <- icell %>% filter(Cell==sp) %>% slice_max(order_by=sagindex,n=100) %>% mutate(sag_group="High")
	icell1 <- icell %>% filter(Cell==sp & sagindex > median(icell_w6$sagindex)) %>% mutate(sag_group="High") 
        icell1 %>% summary()

        icell2 <- icell %>% filter(sagindex >= median100_min & sagindex <= median100_max ) %>% mutate(sag_group="Low")
        icell2 %>% summary()
	
	icell <- rbind(icell1,icell2) %>% as.data.frame()
	tsvn1 <- str_c(pre,"_cell-sag_W3.csv")
        write_tsv(icell,tsvn1)

	# plot
	p <- ggplot(data=icell,aes(sag_group,sagindex,fill=sag_group)) + geom_boxplot() +  theme_cowplot()+  theme(legend.position="none",plot.title = element_text(size=11)) + ggtitle("boxplot of sag-index value for W3  senescence cells")
	ggsave(str_c(pre,"-W3-","boxplot.png"),p, dpi=600)
		
	
	# step 2: prepare rds
	irds <- readRDS(irds)
	rds <- subset(irds,subset=(Cell_ID %in% unique(icell$Cell_ID )))
	
	# step 3: DEGs detection
	tmp <- rds@meta.data
  
  tmp <- tmp %>% left_join(icell,by=c("Cell_ID"="Cell_ID")) %>% select(sag_group)
	
	rds@meta.data$sag_group <- tmp$sag_group
	
	Idents(rds) <- rds@meta.data$sag_group
	
	markers <- FindAllMarkers(rds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	
	markers <- markers %>% group_by(cluster) %>% as.data.frame()
	
	tsvn <- str_c(pre,"_markers_W3.csv")
	write_tsv(markers,tsvn)

	

}

sag_degs34 <-function(args){

	irdsn <- args[1]
	irds <- readRDS(irdsn)

	sags_allf <- args[2]
	#sags_allf <- "/hwfssz1/ST_EARTH/P20Z10200N0035/USER/tancong/project/03.Aradibopsis_scRNA/1.leaf_development/10.leaf_merge/renew_0426/20230501SAGs_index.csv"
	
	outs <- list()
	plots <- list()
	cellTypes<- c("Epidermis","Mesophyll","Vascular")

	for (cellType in cellTypes ){

		# step 1: extract top 10% of cells with sag-index value in stage 4
		
		sags_all <- read_csv(sags_allf) %>% filter(annotation==cellType )
		
		cut_value <- sags_all %>% filter(stage=="Stage_4_rosette") %>% select(sagindex) %>% unlist() %>% quantile(.,0.95)

		sags_all3 <- sags_all %>% filter(stage=="Stage_3_rosette") 

		cells3_high <- sags_all %>% filter(stage=="Stage_3_rosette" & sagindex > cut_value)  %>% mutate(sag_group="High")
		cells3_normal <- sags_all %>% filter(stage=="Stage_3_rosette" & sagindex > quantile(sags_all3$sagindex,0.45) & sagindex < quantile(sags_all3$sagindex,0.55)) %>% mutate(sag_group="Normal")
		
		icell <- rbind(cells3_high,cells3_normal) %>% as.data.frame()
		#cell3_selected <- str_c(pre,"_stage3_selected.csv")
		sheet1n <- str_c(cellType,"S3-selected-cells",sep="_")
		outs[[sheet1n]] <- icell
		#write_tsv(icell,cell3_selected)

		# plot
		plots[[cellType]] <- ggplot(data=icell,aes(sag_group,sagindex,fill=sag_group)) + geom_boxplot(outlier.shape = NA) + ylim(c(0,1.5))+  theme_classic() +  theme(legend.position="none",plot.title = element_text(size=11)) + labs(title=str_c(cellType," SAG "))
		
		# step 2: prepare rds
		
		rds <- subset(irds,cells=icell$cell)
		
		# step 3: DEGs detection
		tmp <- rds@meta.data

		tmp$Cell_ID <- Cells(rds)

	
		tmp <- tmp %>% left_join(icell,by=c("Cell_ID"="cell")) %>% select(sag_group)
		
		rds@meta.data$sag_group <- tmp$sag_group
		
		Idents(rds) <- rds@meta.data$sag_group
		
		markers <- FindAllMarkers(rds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
		
		markers <- markers %>% group_by(cluster) %>% as.data.frame()
		
		sheet2n <- str_c(cellType,"S3-DEGs",sep="_")
		outs[[sheet2n]] <- markers

		#tsvn <- str_c(pre,"_markers_W3.csv")
		#write_tsv(markers,tsvn)

	}


pre <- Sys.Date() %>% str_replace_all("-","") %>% str_c("_S3-comparison_",sep="")

xlsxn <- str_c(pre,".xlsx")
write.xlsx(outs,xlsxn,colnames=T)

pdf(str_c(pre,".pdf"))
print(plots[[cellTypes[1]]] / plots[[cellTypes[2]]] / plots[[cellTypes[3]]])
dev.off()

}


sag_degs44 <-function(args){

	irdsn <- args[1]
	irds <- readRDS(irdsn)

	sags_allf <- args[2]
	#sags_allf <- "/hwfssz1/ST_EARTH/P20Z10200N0035/USER/tancong/project/03.Aradibopsis_scRNA/1.leaf_development/10.leaf_merge/renew_0426/20230501SAGs_index.csv"
	
	outs <- list()
	plots <- list()
	cellTypes<- c("Epidermis","Mesophyll","Vascular")

	for (cellType in cellTypes ){

		# step 1: extract top 10% of cells with sag-index value in stage 4
		
		sags_all <- read_csv(sags_allf) %>% filter(annotation==cellType )
		
		cut_value <- sags_all %>% filter(stage=="Stage_4_rosette") %>% select(sagindex) %>% unlist() %>% quantile(.,0.95)

		sags_all4 <- sags_all %>% filter(stage=="Stage_4_rosette")

		cells4_high <- sags_all %>% filter(stage=="Stage_4_rosette" & sagindex > cut_value)  %>% mutate(sag_group="High")
		cells4_normal <- sags_all %>% filter(stage=="Stage_4_rosette" & sagindex > quantile(sags_all4$sagindex,0.45) & sagindex < quantile(sags_all4$sagindex,0.55)) %>% mutate(sag_group="Normal")
		
		icell <- rbind(cells4_high,cells4_normal) %>% as.data.frame()
		#cell3_selected <- str_c(pre,"_stage3_selected.csv")
		sheet1n <- str_c(cellType,"S4-selected-cells",sep="_")
		outs[[sheet1n]] <- icell
		#write_tsv(icell,cell3_selected)

		# plot
		plots[[cellType]] <- ggplot(data=icell,aes(sag_group,sagindex,fill=sag_group)) + geom_boxplot(outlier.shape = NA) + ylim(c(0,1.5))+  theme_classic() +  theme(legend.position="none",plot.title = element_text(size=11)) + labs(title=str_c(cellType," SAG "))
		
		# step 2: prepare rds
		
		rds <- subset(irds,cells=icell$cell)
		
		# step 3: DEGs detection
		tmp <- rds@meta.data

		tmp$Cell_ID <- Cells(rds)

	
		tmp <- tmp %>% left_join(icell,by=c("Cell_ID"="cell")) %>% select(sag_group)
		
		rds@meta.data$sag_group <- tmp$sag_group
		
		Idents(rds) <- rds@meta.data$sag_group
		
		markers <- FindAllMarkers(rds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
		
		markers <- markers %>% group_by(cluster) %>% as.data.frame()
		
		sheet2n <- str_c(cellType,"S4-DEGs",sep="_")
		outs[[sheet2n]] <- markers

		#tsvn <- str_c(pre,"_markers_W3.csv")
		#write_tsv(markers,tsvn)

	}


pre <- Sys.Date() %>% str_replace_all("-","") %>% str_c("_S4-comparison_",sep="")

xlsxn <- str_c(pre,".xlsx")
write.xlsx(outs,xlsxn,colnames=T)

pdf(str_c(pre,".pdf"))
print(plots[[cellTypes[1]]] / plots[[cellTypes[2]]] / plots[[cellTypes[3]]])
dev.off()

}



sag_degs_s6_late <-function(args){

	
	irdsn <- args[1]
	irds <- readRDS(irdsn)
	sags_allf <- args[2]
	#sags_allf <- "/hwfssz1/ST_EARTH/P20Z10200N0035/USER/tancong/project/03.Aradibopsis_scRNA/1.leaf_development/10.leaf_merge/renew_0426/20230501SAGs_index.csv"
	outs <- list()
	plots <- list()
	cellTypes <-c("Epidermis","Mesophyll","Vascular")

	for (cellType in cellTypes){
		pre <- Sys.Date() %>% str_replace_all("-","") %>% str_c("_sagvalue-comparison",sep="")

		# step 1: extract top 10% of cells with sag-index value in stage 4
		
		sags_all <- read_csv(sags_allf) %>% filter(annotation==cellType & stage=="Stage_6_rosette")

		cells6_sag <- sags_all %>% filter(stage=="Stage_6_rosette") %>% select(sagindex)%>% unlist()
		
		cut_value <- sags_all %>% filter(stage=="Stage_6_rosette") %>% select(sagindex) %>% unlist() %>% quantile(.,0.95)
		#cut_value <- 2.0

		cells6_high <- sags_all %>% filter(stage=="Stage_6_rosette" & sagindex > cut_value)  %>% mutate(sag_group="High")
		cells6_normal <- sags_all %>% filter(stage=="Stage_6_rosette" & sagindex > quantile(cells6_sag,0.45) & sagindex < quantile(cells6_sag,0.55)) %>% mutate(sag_group="Normal")
		# cells3_normal <- sags_all %>% filter(stage=="Stage_3_rosette" & sagindex > quantile(sags_all$sagindex,0.45) & sagindex < quantile(sags_all$sagindex,0.55)) %>% mutate(sag_group="Normal")
		
		icell <- rbind(cells6_high,cells6_normal) %>% as.data.frame()
		sheet1n <- str_c(cellType,"selected-cells",sep="_")
		outs[[sheet1n]] <- icell
		# cell6_selected <- str_c(pre,"_stage6_selected.csv")
		# write_tsv(icell,cell6_selected)

		# plot
		plots[[cellType]] <- ggplot(data=icell,aes(sag_group,sagindex,fill=sag_group)) + geom_boxplot() + ylim(c(0,2.5))+  theme_classic() +  theme(legend.position="none",plot.title = element_text(size=11)) + labs(title=str_c(cellType," SAG "))
		
		#pdf(str_c(pre,"-S6-Sag-compare","boxplot.pdf"))
		#print(p)
		#dev.off()


		# step 2: prepare rds
		#irds <- readRDS(irds)
		rds <- subset(irds,cells=icell$cell)
		
		# step 3: DEGs detection
		tmp <- rds@meta.data

		tmp$Cell_ID <- Cells(rds)

	
		tmp <- tmp %>% left_join(icell,by=c("Cell_ID"="cell")) %>% select(sag_group)
		
		rds@meta.data$sag_group <- tmp$sag_group
		
		Idents(rds) <- rds@meta.data$sag_group
		
		markers <- FindAllMarkers(rds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
		
		markers <- markers %>% group_by(cluster) %>% as.data.frame()
		
		# tsvn <- str_c(pre,"_markers_S6.csv")
		# write_tsv(markers,tsvn)

		sheet2n <- str_c(cellType,"S6-DEGs",sep="_")
		outs[[sheet2n]] <- markers
	}


pre <- Sys.Date() %>% str_replace_all("-","") %>% str_c("_S6-comparison",sep="")

xlsxn <- str_c(pre,".xlsx")
write.xlsx(outs,xlsxn,colnames=T)

pdf(str_c(pre,".pdf"))
print(plots[[cellTypes[1]]] / plots[[cellTypes[2]]] / plots[[cellTypes[3]]])
dev.off()


}

setwd("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/tancong/project/03.Aradibopsis_scRNA/1.leaf_development/10.leaf_merge/renew_0426/")



#irds <- "rosette_combine.rds"
irds <- "/hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhaocaiyao/scRNA_seq/new.data.v1/rossete_all_stage/try3_all_S4_2/harmony_try1_sample1_gene500_8000_lambda1/data.cc.RDS"
sags_allf <- "20230501SAGs_index.csv"

# cellType <- c("Epidermis","Mesophyll","Vascular")
#sag_degs34(c(irds,sags_allf,cellType[1]))
#sag_degs34(c(irds,sags_allf,cellType[2]))


sag_degs34(c(irds,sags_allf))
sag_degs44(c(irds,sags_allf))
#sag_degs_s6_late(c(irds,sags_allf))



