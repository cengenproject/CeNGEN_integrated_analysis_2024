# Scripts and data files for analyses performed in Barrett et al., 2024 from the CeNGEN project

### Scripts are split by figure, and the order to run the scripts is indicated in this readme file

## required packages to install
[LittleBites](https://github.com/alecbarrett/LittleBites)
[Prop2Count](https://github.com/alecbarrett/Prop2Count)

## required files not included here:

"032224_L4_all_cells_Seurat5.rds" found on the [CeNGEN downloads page](https://www.cengen.org/downloads/), the full dataset (100k cells) was used for this analysis. [direct download link](https://cengen.org/storage/032224_L4_all_cells_Seurat5.rds)

### Pre-processing

#### Identifying likely Artifact genes

Some genes were removed from analysis because they are prone to artifactual expression patterns due to their potential inclusion in promoter-fusion constructs used for FACS sorting.

<ins>Putative artifact gene identification described in:</ins> R/filter_marker_promoters.R

### Figure 1

<ins>FACS plots and micrographs are not included here</ins>   
<ins>PCA:</ins> R/plot_PCA.R

### Supplementary Figure 1

<ins>Single Cell - Bulk Correlation plot:</ins> R/Bulk_single_cell_correlation_heatmap.R

### Figure 2

<ins>pre-processing LittleBites:</ins> R/subtraction.R  
<ins>pre-processing bMIND:</ins> R/bMIND_deconvolution.R  
<ins>pre-processing ENIGMA:</ins> R/ ENIGMA_L2_norm.R

<ins>ROC for non-neuronal and ubiquitous genes:</ins> R/UN_ROC_comparison_240401.R  
<ins>ROC for neuronal genes:</ins> R/Neuronal_ROC_comparison_240516.R

### Supplementary Figure 2

<ins>AUROC barcharts for non-neuronal and ubiquitous genes:</ins> R/UN_ROC_comparison_240401.R   
<ins>AUCROC barcharts for neuronal genes:</ins> R/Neuronal_ROC_comparison_240516.R  
<ins>AUROC scatterplots for non-neuronal and ubiquitous genes per sample:</ins> R/UN_ROC_comparison_240401.R  

### Figure 3

<ins>AFD ttx-1 expression:</ins> R/AFD_UMAP.R  
<ins>Single Cell ROC plot:</ins> R/Single_cell_ROC.R  
<ins>Simulated and real logit models:</ins> R/240520_Logit_transform_parameter_comparison.R  

### Supplementary Figure 3

<ins>Running Differential Expression analysis:</ins> R/single_cell_differential_expression.R  
<ins>Visualizing Differential Expression vs Ground Truth:</ins> R/Differential_expression_metrics.R  

### Figure 4

<ins>Generating Integrated Datasets:</ins> R/Integrated_bulk_SingleCell.R

<ins>Neuronal ROC, AUROC, Precision-Recall, and Sensitivity bootstrapping plots:</ins> R/Integrated_Neuronal_ROC.R  
<ins>Barchart of genes detected in single cell and integrated:</ins> R/Integrated_genes_new_detection.R

<ins>GO term and anatomical enrichment performed using the Enrichment Analysis tool on wormbase</ins> https://wormbase.org/tools/enrichment/tea/tea.cgi

### Supplementary Figure 4


### Figure 5

<ins>Scatterplots and bar plots of genes detected in integrated data but missed in single cell:</ins> R/Integrated_genes_new_detection.R

### Supplementary Figure 5

<ins>GO term and anatomical enrichment performed using the Enrichment Analysis tool on wormbase</ins> https://wormbase.org/tools/enrichment/tea/tea.cgi

<ins>Gene lists for each cell that was used as input for the GSEA tool can be extracted using:</ins> R/Integrated_genes_new_detection.R

### Figure 6

<ins>Noncoding RNA pie charts, and specific gene expression pattern heatmap:</ins> R/ncRNA_graphs_and_specificity.R

<ins>GPCR pseudogene expression plots:</ins> R/GPCR_gene_analysis_240731.R

### Supplementary Figure 6

<ins>Relative maximum gene expression distributions, Contamination correlation cutoffs, ubiquitous ncRNA gene distribution, and total ncRNA genes per cell (by gene biotype):</ins> R/ncRNA_graphs_and_specificity.R
