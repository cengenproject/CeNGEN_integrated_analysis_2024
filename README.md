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

Putative artifact gene identification described in: R/filter_marker_promoters.R

### Figure 1

FACS plots and micrographs are not included here
PCA: R/plot_PCA.R

### Supplementary Figure 1

Single Cell - Bulk Correlation plot: R/Bulk_single_cell_correlation_heatmap.R

### Figure 2

pre-processing:
    bMIND: R/bMIND_deconvolution.R
    ENIGMA: R/ ENIGMA_L2_norm.R

ROC for non-neuronal and ubiquitous genes: R/UN_ROC_comparison_240401.R
ROC for neuronal genes: R/Neuronal_ROC_comparison_240516.R

### Supplementary Figure 2

AUROC barcharts for non-neuronal and ubiquitous genes: R/UN_ROC_comparison_240401.R
AUCROC barcharts for neuronal genes: R/Neuronal_ROC_comparison_240516.R
AUROC scatterplots for non-neuronal and ubiquitous genes per sample: R/UN_ROC_comparison_240401.R

### Figure 3

AFD ttx-1 expression: R/AFD_UMAP.R
Single Cell ROC plot: R/
Simulated and real logit models: R/240520_Logit_transform_parameter_comparison.R

### Supplementary Figure 3

Running Differential Expression analysis: R/single_cell_differential_expression.R
Visualizing Differential Expression vs Ground Truth: R/Differential_expression_metrics.R

### Figure 4

Generating Integrated Datasets: R/Integrated_bulk_SingleCell.R

Neuronal ROC, AUROC, Precision-Recall, and Sensitivity bootstrapping plots: R/Integrated_Neuronal_ROC.R
Barchart of genes detected in single cell and integrated: R/Integrated_genes_new_detection.R

GO term and anatomical enrichment performed using the Enrichment Analysis tool on wormbase https://wormbase.org/tools/enrichment/tea/tea.cgi

### Supplementary Figure 4


### Figure 5

Scatterplots and bar plots of genes detected in integrated data but missed in single cell: R/Integrated_genes_new_detection.R

### Supplementary Figure 5

GO term and anatomical enrichment performed using the Enrichment Analysis tool on wormbase https://wormbase.org/tools/enrichment/tea/tea.cgi

Gene lists for each cell that was used as input for the GSEA tool can be extracted using: R/Integrated_genes_new_detection.R

### Figure 6

Noncoding RNA pie charts, and specific gene expression pattern heatmap: R/ncRNA_graphs_and_specificity.R

GPCR pseudogene expression plots: R/GPCR_gene_analysis_240731.R

### Supplementary Figure 6

Relative maximum gene expression distributions, Contamination correlation cutoffs, ubiquitous ncRNA gene distribution, and total ncRNA genes per cell (by gene biotype): R/ncRNA_graphs_and_specificity.R
