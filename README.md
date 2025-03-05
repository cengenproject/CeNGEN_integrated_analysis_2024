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

<img width="300" alt="Figure 1B" src="https://github.com/user-attachments/assets/bda05119-cf47-412a-9dcd-fcbeae0162b0" />

### Supplementary Figure 1

<ins>Single Cell - Bulk Correlation plot:</ins> R/Bulk_single_cell_correlation_heatmap.R

<img width="300" alt="Figure S1B" src="https://github.com/user-attachments/assets/273afe18-7df7-4519-9165-60fdfd528985" />

### Figure 2

<ins>pre-processing LittleBites:</ins> R/subtraction.R  
<ins>pre-processing bMIND:</ins> R/bMIND_deconvolution.R  
<ins>pre-processing ENIGMA:</ins> R/ ENIGMA_L2_norm.R

<ins>ROC for non-neuronal and ubiquitous genes:</ins> R/UN_ROC_comparison_240401.R  

<img width="394" alt="Figure 2B" src="https://github.com/user-attachments/assets/410afa91-0f1d-4542-b0bc-147dd81d1dfe" />

<ins>ROC for neuronal genes:</ins> R/Neuronal_ROC_comparison_240516.R

<img width="396" alt="Figure 2C" src="https://github.com/user-attachments/assets/66c44c39-a7b0-4e3b-8be8-400c355525dc" />

### Supplementary Figure 2

<ins>AUROC barcharts for non-neuronal and ubiquitous genes:</ins> R/UN_ROC_comparison_240401.R   

<img width="398" alt="Figure S2A" src="https://github.com/user-attachments/assets/82e163e4-3438-4b81-bf24-6950e3f39bd7" />

<ins>AUCROC barcharts for neuronal genes:</ins> R/Neuronal_ROC_comparison_240516.R  

<img width="300" alt="Figure S2B" src="https://github.com/user-attachments/assets/f351dfd3-cc28-45f3-ab68-854046ed03e4" />

<ins>AUROC scatterplots for non-neuronal and ubiquitous genes per sample:</ins> R/UN_ROC_comparison_240401.R  

<img width="300" alt="" src="https://github.com/user-attachments/assets/78d10ce4-8afa-4167-a77a-d99e703f73d9" />

### Figure 3

<ins>AFD ttx-1 expression:</ins> R/AFD_UMAP.R  

<img width="300" alt="" src="https://github.com/user-attachments/assets/8a91eac7-93d8-4573-a99e-c8c4839bee6c" />

<ins>Single Cell ROC plot:</ins> R/Single_cell_ROC.R

<img width="300" alt="" src="https://github.com/user-attachments/assets/0c0a5b26-f7fd-4f13-81ed-2f9e19286444" />

<ins>Simulated and real logit models:</ins> R/240520_Logit_transform_parameter_comparison.R  

<img width="300" alt="" src="https://github.com/user-attachments/assets/7ccbce48-2c63-45a2-a43c-55b240892deb" />

### Supplementary Figure 3

<ins>Running Differential Expression analysis:</ins> R/single_cell_differential_expression.R  
<ins>Visualizing Differential Expression vs Ground Truth:</ins> R/Differential_expression_metrics.R  
<img width="300" alt="Screenshot 2025-03-04 at 9 16 27 PM" src="https://github.com/user-attachments/assets/e290d2d3-6b3a-40f1-bc01-8d0e1f4ea2bd" />

### Figure 4

<ins>Generating Integrated Datasets:</ins> R/Integrated_bulk_SingleCell.R

<ins>Neuronal ROC, AUROC, Precision-Recall, and Sensitivity bootstrapping plots:</ins> R/Integrated_Neuronal_ROC.R  
<ins>Barchart of genes detected in single cell and integrated:</ins> R/Integrated_genes_new_detection.R

<img width="400" alt="" src="https://github.com/user-attachments/assets/5361a251-e795-4e10-ba48-b4bb8e035cac" />

<ins>GO term and anatomical enrichment performed using the Enrichment Analysis tool on wormbase</ins> https://wormbase.org/tools/enrichment/tea/tea.cgi

<img width="300" alt="" src="https://github.com/user-attachments/assets/dfbf42dd-8bca-4c0f-a29b-5b2ca2cbb685" />



### Supplementary Figure 4

<ins>Heatmaps of expresssion vs ground truth correspondence</ins>

<img width="400" alt="" src="https://github.com/user-attachments/assets/ef4c02cf-5ed6-41b8-9cc8-1c05416cec84" />

### Figure 5

<ins>Scatterplots and bar plots of genes detected in integrated data but missed in single cell:</ins> R/Integrated_genes_new_detection.R

<img width="300" alt="" src="https://github.com/user-attachments/assets/2d96637b-cfdc-4b54-8974-2e8071f2d141" />

### Supplementary Figure 5

<ins>GO term and anatomical enrichment performed using the Enrichment Analysis tool on wormbase</ins> https://wormbase.org/tools/enrichment/tea/tea.cgi

<ins>Gene lists for each cell that was used as input for the GSEA tool can be extracted using:</ins> R/Integrated_genes_new_detection.R

<img width="300" alt="" src="https://github.com/user-attachments/assets/8ee4b9ec-0a90-4734-b4ef-d57653cb2156" />

### Figure 6

<ins>Noncoding RNA pie charts, and specific gene expression pattern heatmap:</ins> R/ncRNA_graphs_and_specificity.R

<img width="300" alt="" src="https://github.com/user-attachments/assets/2fc26708-ccc1-43c1-a434-cfaec5a5d4a8" />

<ins>GPCR pseudogene expression plots:</ins> R/GPCR_gene_analysis_240731.R

<img width="300" alt="" src="https://github.com/user-attachments/assets/9e03e36b-ba11-4fb0-965d-af25ab739f80" />

### Supplementary Figure 6

<ins>Relative maximum gene expression distributions, Contamination correlation cutoffs, ubiquitous ncRNA gene distribution, and total ncRNA genes per cell (by gene biotype):</ins> R/ncRNA_graphs_and_specificity.R

<img width="300" alt="" src="https://github.com/user-attachments/assets/31038fc6-42d5-49b8-8042-866152727f9d" />
<img width="300" alt="" src="https://github.com/user-attachments/assets/c593e6c3-1754-4aed-bff1-069502e9fcd5" />

