
#### running through ENIGMA testing

# install.packages(c("Matrix","S4Vectors","corpcor","MASS","e1071","ggplot2","cowplot","magrittr","purrr","tibble","nnls","doParallel","tidyr","plyr","vctrs","matrixStats"))
# BiocManager::install(c("SingleCellExperiment","scater","Biobase","SummarizedExperiment","sva","preprocessCore"))
# devtools::install_github("WWXKenmo/ENIGMA_test")


library(ggplot2)
library(ENIGMA)
library(stringr)
library(ComplexHeatmap)


## load in bulk data

Data <- read.table('Data/Bulk_data_bsn12_231211.tsv')

Data <- Data[apply(Data, 1, sd) > 0,] ### remove invariant genes

cell_types <- unique(str_split_fixed(colnames(Data), 'r', 2)[,1]) ## get a list of cell types


## load in single cell prop2count data

seurat_pseudobulk_tissue <- read.table('Data/Tissue_pseudobulk_prop2count_082022.tsv')


sc_meta <- data.frame(row.names = colnames(seurat_pseudobulk_tissue),
                      sample = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,2],
                      cell_type = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,1])


seurat_pseudobulk_tissue_CPM <-  sweep(seurat_pseudobulk_tissue, MARGIN = 2, STATS = colSums(seurat_pseudobulk_tissue), FUN = '/') * 1000000

seurat_pseudobulk_tissue_ave_CPM <- sapply(unique(sc_meta$cell_type), function(x){
  tmp <- seurat_pseudobulk_tissue_CPM[, str_split_fixed(colnames(seurat_pseudobulk_tissue_CPM), '__', 2)[,1]==x]
  return(rowMeans(tmp))
})



#### subset to common genes

enigma.common.genes <- intersect(rownames(Data), rownames(seurat_pseudobulk_tissue))

Data <- Data[enigma.common.genes,]


#### make the ENIGMA object
egm_aggre <- ENIGMA::create_ENIGMA(bulk = as.matrix(Data), ref = as.matrix(seurat_pseudobulk_tissue_ave_CPM), ref_type = "aggre")

#egm_aggre <- batch_correct(egm_aggre) ## batch correction gives poor AUROC compared to no batch correction

egm_aggre <- get_cell_proportion(egm_aggre, method = "RLR")

Heatmap(egm_aggre@result_cell_proportion |> t(), cluster_columns = F)

#egmL2 <- ENIGMA_L2_max_norm(egm_aggre, preprocess = 'sqrt') ## sqrt preprocessing does not perform as well as log2(+1) on gene detection
egmL2_log <- ENIGMA_L2_max_norm(egm_aggre, preprocess = 'log')

saveRDS(egmL2_log, 'Data_out/bulk_bsn12_enigma_L2_norm_log_deconv_030824.rds')


## extract the neuron specific columns
egmL2_log_neuron <- egmL2_log@result_CSE@assays@data$counts[, egmL2_log@result_CSE$cell_type=='Neuron']
colnames(egmL2_log_neuron) <- str_split_fixed(colnames(egmL2_log_neuron), ':', 2)[,1]


write.table(egmL2_log_neuron, 'Data_out/bulk_bsn12_enigma_L2_norm_log_NeuronalCounts_030824.tsv', sep = '\t')

