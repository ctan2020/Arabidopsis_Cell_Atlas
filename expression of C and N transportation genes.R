
#

library(tidyverse)
library(ggplotify)
library(openxlsx)
library(patchwork)
library(cowplot)
library(ggpubr)
library(readxl)
library(pheatmap)


setwd('../path/work-directory/')

ifile1 <- "rosette_exp_0710.xlsx"

leaf_gene_exp <- read_xlsx(ifile1,col_names = T) %>% as.data.frame()  %>% select(!contains(c("base","tip","mid")))  %>%
  select("Gene_ID",contains(c("Epi","Mes","PP","CC")))

leaf_gene_exp %>% names() 

ifile2 <- "flower_silique_selected.xlsx"

flower_gene_exp <- read_xlsx(ifile2,col_names = T) %>% as.data.frame() %>% select(c(1,3,2,6,5,4,8,7,11,10,9))

flower_gene_exp %>% names() 

ifile3 <- "CN-gene-list-0823-1802.xlsx"  
  
C_genes_all <- read_xlsx(ifile3,col_names = TRUE,sheet = "C-genes-all") %>% as.data.frame() 
C_genes_all$Gene.Locus <- C_genes_all$Gene.Locus %>% str_to_upper()

C_genes_ordered <-  read_xlsx(ifile3,col_names = TRUE,sheet = "C-genes-ordered") %>% as.data.frame() 
C_genes_ordered$Gene.Locus <- C_genes_ordered$Gene.Locus %>% str_to_upper()


C_genes_all_leaf_exp <- left_join(C_genes_all,leaf_gene_exp,by=c("Gene.Locus"="Gene_ID"))
C_genes_all_flower_exp <- left_join(C_genes_all,flower_gene_exp,by=c("Gene.Locus"="Gene.Locus"))

exp <-C_genes_all_leaf_exp

anno_col <- data.frame(Cell.Type=factor(rep(c("EPI","MES","PP","CC"),each=6),levels = c("EPI","MES","PP","CC")))
row.names(anno_col) <- colnames(select(exp,-seq(4)))

p2 <- pheatmap(select(exp,-seq(4)),scale = "row",cluster_rows = F,cluster_cols = F,
               cellheight = 6,
               labels_row = str_c(exp$Gene.Locus," | ",exp$Gene.Name),
               gaps_col = c(6,12,18),fontsize_row = 6,border_color = "gray",
               annotation_col = anno_col,labels_col =rep(str_c(rep("S",6),seq(1,6)),4),
               angle_col = 45)

g2 <- as.ggplot(p2)

g2

pre <- Sys.Date() %>% str_replace_all("-","")
ofile <- str_c(pre,"_C_genes_all_leaf_exp.pdf")
ggsave(ofile,g2,width = 10,height = 6)


# heatmap for c pathways of flower
exp <- C_genes_all_flower_exp
anno_col <- data.frame(Stage=factor(rep(c('Flower(early stage)','Silique D5'),each=5)),
                       Cell_Type=factor(c("PP","CC","Nectary","Tapetum","Pollen","PP","CC","Seed coat","Endosperm","Embryo"),levels=c("PP","CC","Nectary","Tapetum","Pollen","Seed coat","Endosperm","Embryo")))
row.names(anno_col) <- colnames(select(exp,-seq(4)))


p <- pheatmap(select(exp,-seq(4)),scale = "row",cluster_rows = F,cluster_cols = F,
               cellheight = 6,
               labels_row = str_c(exp$Gene.Locus," | ",exp$Gene.Name),
               gaps_col = c(5),fontsize_row = 6,border_color = "gray",
              annotation_col = anno_col,labels_col = anno_col$Cell_Type,
              angle_col = 45)

g <- as.ggplot(p)

g

pre <- Sys.Date() %>% str_replace_all("-","")
ofile <- str_c(pre,"_C_genes_all_flower_exp.pdf")
ggsave(ofile,g,width = 10,height = 8)

N_genes_all<-  read_xlsx(ifile3,col_names = TRUE,sheet = "N-genes-all") %>% as.data.frame() 
N_genes_ordered <-  read_xlsx(ifile3,col_names = TRUE,sheet = "N-genes-ordered") %>% as.data.frame() 
N_genes_ordered$Gene.Locus <- N_genes_ordered$Gene.Locus %>% str_to_upper()
N_genes_all$Gene.Locus <- N_genes_all$Gene.Locus %>% str_to_upper()

N_genes_ordered_leaf_exp <- left_join(N_genes_ordered,leaf_gene_exp,by=c("Gene.Locus"="Gene_ID"))
N_genes_ordered_flower_exp <- left_join(N_genes_ordered,flower_gene_exp,by=c("Gene.Locus"="Gene.Locus"))
N_genes_all_leaf_exp <- left_join(N_genes_all,leaf_gene_exp,by=c("Gene.Locus"="Gene_ID"))
N_genes_all_flower_exp <- left_join(N_genes_all,flower_gene_exp,by=c("Gene.Locus"="Gene.Locus"))

# heatmap for 

# heatmap for flower N ordered genes
exp <- N_genes_all_flower_exp

anno_col <- data.frame(Stage=factor(rep(c('Flower(early stage)','Silique D5'),each=5)),
                       Cell_Type=factor(c("PP","CC","Nectary","Tapetum","Pollen","PP","CC","Seed coat","Endosperm","Embryo"),levels=c("PP","CC","Nectary","Tapetum","Pollen","Seed coat","Endosperm","Embryo")))
row.names(anno_col) <- colnames(select(exp,-seq(4)))

p <- pheatmap(select(exp,-seq(4)),scale = "row",cluster_rows = F,cluster_cols = F, 
          cellheight = 6,
          labels_row = str_c(exp$Gene.Locus," | ",exp$Gene.Name),
          gaps_col = c(5),fontsize_row = 6,border_color = "gray",
          annotation_col = anno_col,labels_col = anno_col$Cell_Type,
          angle_col = 45)     
g1 <- as.ggplot(p)
g1


pre <- Sys.Date() %>% str_replace_all("-","")
ofile <- str_c(pre,"_N_genes_all_flower_exp.pdf")

ggsave(ofile,g1,width = 10,height = 18)

#dev.off()

# heatmap for all N genes in  leaf
exp <- N_genes_all_leaf_exp

anno_col <- data.frame(Cell.Type=factor(rep(c("EPI","MES","PP","CC"),each=6),levels = c("EPI","MES","PP","CC")))
row.names(anno_col) <- colnames(select(exp,-seq(4)))


#pdf(paste0("C_gene_ordered_0730","_heatmap.pdf"),height = 5,width = 10)
p2 <- pheatmap(select(exp,-seq(4)),scale = "row",cluster_rows = F,cluster_cols = F,
               cellheight = 6,
               labels_row = str_c(exp$Gene.Locus," | ",exp$Gene.Name),
               gaps_col = c(6,12,18),fontsize_row = 6,border_color = "gray",
               annotation_col = anno_col,labels_col =rep(str_c(rep("S",6),seq(1,6)),4),
               angle_col = 45)

g2 <- as.ggplot(p2)

g2

pre <- Sys.Date() %>% str_replace_all("-","")
ofile <- str_c(pre,"_N_genes_all_leaf_exp.pdf")

ggsave(ofile,g2,width = 10,height = 15)

out <- list()

out[["C_pathway_leaf_exp"]] <- C_genes_all_leaf_exp[which(C_genes_all_leaf_exp$Epidermis_S1>0),]
out[["C_pathway_flower_silique_exp"]] <- C_genes_all_flower_exp[which(C_genes_all_flower_exp$early_flower_phloem.parenchyma>0),]

out[["N_pathway_leaf_exp"]] <- N_genes_all_leaf_exp[which(N_genes_all_leaf_exp$Epidermis_S1>0),]
out[["N_pathway_flower_silique_exp"]] <- N_genes_all_flower_exp[N_genes_all_flower_exp$early_flower_phloem >0 ,]

ofile <- Sys.Date() %>% str_replace_all("-","") %>% str_c("_C-N-pathway-genes-leaf-flower-silique.xlsx")
write.xlsx(out,ofile )

##################################################################







### for ordered genes


exp <- N_genes_ordered_leaf_exp

anno_col <- data.frame(Cell.Type=factor(rep(c("EPI","MES","PP","CC"),each=6),levels = c("EPI","MES","PP","CC")))
row.names(anno_col) <- colnames(select(exp2,-seq(2)))


#pdf(paste0("C_gene_ordered_0730","_heatmap.pdf"),height = 5,width = 10)

p2 <- pheatmap(select(exp2,-seq(2)),scale = "row",cluster_rows = F,cluster_cols = F,
         cellheight = 10,cellwidth = 20,labels_row = str_c(exp2$Gene.Locus," | ",exp2$Gene.Name),
         gaps_col = c(6,12,18),annotation_col = anno_col,labels_col =rep(str_c(rep("S",6),seq(1,6)),4),angle_col = 45)

g2 <- as.ggplot(p2)

ggsave("test.pdf",g2+g1,width = 20,height = 10)

################### for c genes


C_genes_ordered_leaf_exp <- left_join(C_genes_ordered,leaf_gene_exp,by=c("Gene.Locus"="Gene_ID"))
C_genes_ordered_flower_exp <- left_join(C_genes_ordered,flower_gene_exp,by=c("Gene.Locus"="Gene.Locus"))


# heatmap for flower 
exp <- C_genes_ordered_flower_exp
anno_col <- data.frame(Stage=factor(rep(c('Flower(early stage)','Silique D5'),each=5)),
                       Cell_Type=factor(c("PP","CC","Nectary","Tapetum","Pollen","PP","CC","Seed coat","Endosperm","Embryo"),levels=c("PP","CC","Nectary","Tapetum","Pollen","Seed coat","Endosperm","Embryo")))
row.names(anno_col) <- colnames(select(exp,-seq(2)))


#pdf(paste0("leaf_Ngenes_selected0731","_heatmap.pdf"),height = 25,width = 20)
#pheatmap(exp2,scale = "row",cluster_rows = F,cluster_cols = F,labels_row = unique(gene_groups$Labels[which(gene_groups$Gene.Locus %in% exp$Gene_ID)]),gaps_col = c(6,12,18) )
#pheatmap(select(exp,-1),scale = "row",cluster_rows = F,cluster_cols = F ,
#        gaps_col = c(6,12,18),gaps_row=seq(0,150,by=20))

p1 <- pheatmap(select(exp,-seq(2)),scale = "row",cluster_rows = F,cluster_cols = F,
               gaps_col = c(5),cellwidth = 20,
               annotation_col = anno_col,labels_col = anno_col$Cell_Type,
               #     labels_col =rep(str_c(rep("S",6),seq(1,6)),4),
               angle_col = 45,labels_row = str_c(exp$Gene.Locus," | ",exp$Gene.Name),cellheight = 10)

g1 <- as.ggplot(p1)

#dev.off()

# heatmap for C gene of leaf
exp2 <- C_genes_ordered_leaf_exp

anno_col <- data.frame(Cell.Type=factor(rep(c("EPI","MES","PP","CC"),each=6),levels = c("EPI","MES","PP","CC")))
row.names(anno_col) <- colnames(select(exp2,-seq(2)))


#pdf(paste0("C_gene_ordered_0730","_heatmap.pdf"),height = 5,width = 10)

p2 <- pheatmap(select(exp2,-seq(2)),scale = "row",cluster_rows = F,cluster_cols = F,
               cellheight = 10,cellwidth = 20,labels_row = str_c(exp2$Gene.Locus," | ",exp2$Gene.Name),
               gaps_col = c(6,12,18),annotation_col = anno_col,labels_col =rep(str_c(rep("S",6),seq(1,6)),4),angle_col = 45)

g2 <- as.ggplot(p2)

ggsave("test.pdf",g2+g1,width = 20,height = 10)






#dev.off()


########################



 #exp <- all_gene_exp %>% filter(Gene_ID %in% genes) 
 #exp$Gene_ID <- factor(exp$Gene_ID, levels=gene_groups$Gene.Locus,labels=gene_groups$Labels )
 exp <- left_join(gene_groups,all_gene_exp,by=c("Gene.Locus"="Gene_ID"))
 
  
  anno_col <- data.frame(Cell.Type=factor(rep(c("EPI","MES","PP","CC"),each=6),levels = c("EPI","MES","PP","CC")))
  row.names(anno_col) <- colnames(select(exp,-seq(4)))
  
  
pdf(paste0("C_gene_ordered_0730","_heatmap.pdf"),height = 5,width = 10)
  
pheatmap(select(exp,-seq(4)),scale = "row",cluster_rows = F,cluster_cols = F,cellheight = 10,cellwidth = 20,labels_row = str_c(exp$Gene.Locus," | ",exp$Gene.Name),
         gaps_col = c(6,12,18),annotation_col = anno_col,labels_col =rep(str_c(rep("S",6),seq(1,6)),4),angle_col = 45)

    dev.off()
  
# }
  
###

p1 = as.ggplot(pheatmap(select(exp,-1),scale = "row",cluster_rows = F,cluster_cols = F ,
                         gaps_col = c(6,12,18),gaps_row=seq(0,150,by=20)))
print(p1)
