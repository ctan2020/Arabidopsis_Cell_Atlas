
An Arabidopsis single-nucleus atlas decodes leaf senescence and nutrient allocation

# 1. Summary

Over one million nuclei from 20 samples passed accepted droplet-based single nuclei filtering criteria with gene number > 200. Following the majority of previous studies, we used a more rigorous criteria with gene number > 500 to obtain high-quality nuclei for downstream analyses. A total of 913,769 nuclei passed our quality control, with 1610 mean genes and 2451 mean UMIs captured in each nucleus. Cells from each of the 20 samples were independently clustered using Seurat and integrated to generate a comprehensive transcriptome atlas. Independent clustering of each dataset resulted in a total of 367 clusters, 97.55% (358 out of 367) of which were successfully annotated to corresponding cell types based on known cell-type marker genes or biological functions of cluster-enriched genes.

This atlas serves as a valuable resource for the entire Arabidopsis community and contributes a comprehensive reference for future investigations of non-model plant species at single cell resolution.
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/Fig1-0828-01.png" width="800" high="10400"/></center></td>

Figure 1 | A comprehensive cell atlas of Arabidopsis.
a. Sampling design of Arabidopsis organs analysed in this study. A total of 20 samples at different developmental stages were collected from individuals cross 6 day to 49 day. 
b. Numbers of nuclei profiled and genes for each sample captured. 
c. A UMAP visualization is shown of global clustering of all cells from the dataset colored by organ. 
b.  UMAP visualization of all clusters of global cell atlas colored by major cell types.

# 2. Data description
## 2.1 Raw data
In total, 20 samples were collected at indicated day post-germination for all tissues including seedling, cotyledon, hypocotyl, root, rosette leaf, stem, lateral leaf, flower, and silique. Then we used an in-house nuclei isolation protocol (see method) and DNBelab C Series Single-Cell Library Prep Set (MGI, 1000021082) for snRNA data generation as previously described (Han et al., 2022). The concentration of DNA library was measured by Qubit (Invitrogen). Libraries were sequenced by DNBSEQ-T7RS.
## 2.2 Expression matrix
 The raw sequencing reads were filtered and demultiplexed by PISA (Version 1.1.0) (https://github.com/shiquan/PISA), and aligned to the TAIR10 reference genome using STAR (version 2.7.4a) (Dobin et al., 2013) with default parameters. Over one million nuclei from 20 samples passed the accepted droplet-based single nuclei filtering criteria with > 200 genes/nucleus. Following the majority of previous studies, we used a more rigorous criteria with > 500 genes/nucleus to obtain high-quality nuclei for downstream analyses. A total of 913,769 nuclei passed our quality control, with 1610 mean genes and 2451 mean UMIs captured in each nucleus. Cells from each of the 20 samples were independently clustered to obtain single-sample maps (Extended Data Figs. 1–20), then integrated to generate multi-stage organ level maps and a comprehensive all-cell transcriptome atlas.

# 3. Results

## 3.1 Root (6 DAG)

<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-01.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 1 | SnRNA-seq profiling of root (6 DAG). 

UMAP clustering of root (6 DAG) and annotation of cell types.

## 3.2 Root (11 DAG)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-02.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 2 | SnRNA-seq profiling of root (11 DAG).

UMAP clustering of root (11 DAG) and annotation of cell types.

## 3.3 Shoot (6 DAG)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-03.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 3 | SnRNA-seq profiling of shoot (6 DAG).

UMAP clustering of shoot (6 DAG) and annotation of cell types.

## 3.4 Stem
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-04.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 4 | SnRNA-seq profiling of stem.

UMAP clustering of stem and annotation of cell types.

## 3.5 Cauline
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-05.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 5 | SnRNA-seq profiling of cauline. 

UMAP clustering of cauline and annotation of cell types.

## 3.6 Rosette at stage 1 (14 DAG)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-06.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 6 | SnRNA-seq profiling of rosette at stage 1 (14 DAG).

UMAP clustering of rosette at stage 1 (14 DAG) and annotation of cell types.

## 3.7 Rosette at stage 2 (21DAG)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-07.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 7 | SnRNA-seq profiling of rosette at stage 2 (21 DAG).

UMAP clustering of rosette at stage 2 (21 DAG) and annotation of cell types.

## 3.8 Rosette at stage 3 (28DAG)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-08.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 8 | SnRNA-seq profiling of rosette at stage 3 (28 DAG).

UMAP clustering of rosette at stage 3 (28 DAG) and annotation of cell types.

## 3.9 Rosette at stage 4 (35 DAG)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-09.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 9 | SnRNA-seq profiling of rosette at stage 4 (35 DAG).

UMAP clustering of rosette at stage 4 (35 DAG) and annotation of cell types.

## 3.10 Rosette at stage 5 (42 DAG)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-10.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 10 | SnRNA-seq profiling of rosette at stage 5 (42 DAG).

UMAP clustering of rosette at stage 5 (42 DAG) and annotation of cell types.

## 3.11 Rosette at stage 6 (49 DAG)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-11.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 11 | SnRNA-seq profiling of rosette at stage 6 (49 DAG).

UMAP clustering of rosette at stage 6 (49 DAG) and annotation of cell types.

## 3.12 Early flower
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-12.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 12 | SnRNA-seq profiling of early flower. 

UMAP clustering of early flower and annotation of cell types.

## 3.13 Middle flower
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-13.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 13 | SnRNA-seq profiling of middle flower. 

UMAP clustering of middle flower and annotation of cell types.

## 3.14 Late flower
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-14.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 14 | SnRNA-seq profiling of late flower. 

UMAP clustering of late flower and annotation of cell types.

## 3.15 Silique (0 DPA)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-15.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 15 | SnRNA-seq profiling of silique (0 DPA).

UMAP clustering of silique (0 DPA) and annotation of cell types.

## 3.16 Silique (1 DPA)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-16.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 16 | SnRNA-seq profiling of silique (1 DPA).

UMAP clustering of silique (1 DPA) and annotation of cell types.

## 3.17 Silique (2 DPA)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-17.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 17 | SnRNA-seq profiling of silique (2 DPA).

UMAP clustering of silique (2 DPA) and annotation of cell types.

## 3.18 Silique (3 DPA)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-18.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 18 | SnRNA-seq profiling of silique (3 DPA).

UMAP clustering of silique (3 DPA) and annotation of cell types.

## 3.19 Silique (4 DPA)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-19.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 19 | SnRNA-seq profiling of silique (4 DPA).

UMAP clustering of silique (4 DPA) and annotation of cell types.

## 3.20 Silique (5 DPA)
<td><center><img src="https://southpublic-database.obs.cn-south-1.myhuaweicloud.com/Arabidopsis/figure/20 samples umap figures-20.jpg" width="800" high="10400"/></center></td>
Extended Data Fig. 20 | SnRNA-seq profiling of silique (5 DPA).

UMAP clustering of silique (5 DPA) and annotation of cell types.

