
library(tidyverse)
library(Seurat)
library(getopt)
library(cowplot)
library(patchwork)



arg <- matrix(c("help","h","0","logical","help information","",
                "input","i","1","character","query rds","",
                "input2","f","1","character","query rds","",
                "work_dir","w","1","character","work directory","",
                "sample_id","s","1","character","sample id",""
                ),byrow=T,ncol=6)

opt = getopt(arg)

if(!is.null(opt$help) || is.null(opt$input)){
    cat(paste(getopt(arg, usage = T), "\n"))
    q()
}

setwd(opt$work_dir)

pre <- Sys.Date() %>% str_replace_all("-","") %>% str_c("_",opt$sample_id,"_")



SAG_INDEX <- function(rds1){

        # calculate average expression of given genes at stage 5 and 6 
	
	rds1 <- readRDS(rds1)

        rds1meta <- rds1@meta.data %>% as.data.frame() %>% select(Cell_ID,Stage)

        data1 <- rds1@assays$RNA@data %>% as.data.frame() %>% 
                subset(rownames(.) %in% genes$Gene_ID)
                
        cell56 <- rds1meta %>% filter(Stage =="W5"|Stage =="W6") %>% rownames()

        cell56exp <- data1[,names(data1) %in% cell56]

        geneAve <- apply(cell56exp,1,mean)
        

        # calculate SAG index for each cell

        rela_exp <- apply(data1,2,function(x){return(x/geneAve)}); 

        sag_index <- apply(na.omit(rela_exp),2,mean)
        
        sag_index <- data.frame(sagindex=sag_index)

        tmp <- merge(rds1meta,sag_index,by='row.names')
        
        return(tmp)
}


#write_csv(tmp,"tmp-sagindex.csv")

rdsf <- read_tsv(opt$input)

genes <- read_tsv(opt$input2)

df <- data.frame()
for (one in 1:3 ){
	sag_index <- SAG_INDEX(rdsf$Rds[one])
	sag_index$Cell <- rdsf$Cell[one]
	df <- bind_rows(df,sag_index)	
}

ofile <- str_c(pre,".csv")
write_csv(df,ofile)

##########################

library(tidyverse)
library(forcats)
library(cowplot)

md <- read_csv("20211011_SAG-index.csv")

ggplot(md,aes(Stage,sagindex,fill=Stage)) + 
	geom_boxplot() + ylim(0,10) + 
	geom_jitter(width=0.1,alpha=0.2) +
	facet_grid(~Cell) + theme_cowplot() +   
	ggsave("test.pdf",width=8,height=6)



