
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
