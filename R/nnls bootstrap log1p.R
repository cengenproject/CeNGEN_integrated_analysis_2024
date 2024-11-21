library(Seurat)
library(pbapply)
library(tidyverse)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(nnls)

### seurat object downloaded from CeNGEN website downloads page, ~0.5 gb disk space, ~ 2gb RAM to open
sc_object <- readRDS('~/Bioinformatics/single_cell_data/100720_L4_all_cells_Seurat.rds')

sc_object_counts <- sc_object@assays$RNA@counts
sc_object_meta <- sc_object@meta.data

sc_object <- CreateSeuratObject(counts = sc_object_counts, meta.data = sc_object_meta)

rownames(sc_object@assays$RNA@features@.Data) #<- rownames(sc_object_counts)


sc_object <- sc_object[,sc_object$Tissue != 'Unknown' & sc_object$Tissue != 'Unannotated']


sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('DB01')] <- 'DB'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('DA9')] <- 'DA'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VC_4_5')] <- 'VC'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VB01', 'VB02')] <- 'VB'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('RMD_DV', 'RMD_LR')] <- 'RMD'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('RME_DV', 'RME_LR')] <- 'RME'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VA12')] <- 'VA'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('IL2_DV', 'IL2_LR')] <- 'IL2'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('AWC_ON', 'AWC_OFF')] <- 'AWC'
sc_object <- sc_object[,!(sc_object$Cell.type %in% c('RIV_stressed', 'SMD_stressed'))]


non_neuronal_list <- c('Excretory', 'Glia', 'Hypodermis', 'Intestine', 'Muscle_mesoderm', 'Pharynx', 'Reproductive')


sc_object$neuron_level <- sc_object$Tissue

sc_object$neuron_level[sc_object$neuron_level=='Neuron'] <- sc_object$Cell.type[sc_object$neuron_level=='Neuron']

sc_object_neuron <- sc_object[,sc_object$Tissue=='Neuron']

sc_object_neuron$replicate <- paste(sc_object_neuron$Cell.type, sc_object_neuron$Experiment, sep = '__')

replicate_size <- sapply(unique(sc_object_neuron$replicate), function(replica){
  sum(sc_object_neuron$replicate == replica)
})
replicate_size <- replicate_size[order(names(replicate_size))]
write.table(replicate_size, 'CeNGEN_SC_Biorep_size_061422.tsv', sep = '\t')


sc_size <- sapply(unique(sc_object$neuron_level), function(cell){
  sum(sc_object$neuron_level == cell)
})

sc_size <- sc_size[order(sc_size)]
sc_size['DD'] <- sc_size['VD_DD']
sc_size['VD'] <- sc_size['VD_DD']


write.table(sc_size, 'sc_size_032322.csv', sep = ',')

## subset to 12 cells maximum


### cut down the size of the Seurat object to dramatically speed up the for loop
sc_object_cut <- CreateSeuratObject(counts = sc_object@assays$RNA@layers$counts, 
                                    meta.data = sc_object@meta.data[,c('neuron_level', 'Cell.type', 'orig.ident', 'Barcode', 'Tissue', 'Experiment')])

sc_object_cut@assays$RNA@features@.Data |> rownames() <- sc_object@assays$RNA@features@.Data |> rownames()

sc_object_cut$replicate <- paste(sc_object_cut$neuron_level, sc_object_cut$orig.ident, sep = "__")



NNLS_30_list_log1p <- pblapply(seq(101,200,1), function(seed){
  set.seed(seed)
  to_keep <- lapply(unique(sc_object_cut$neuron_level), function(cell){ ### list of cell barcodes to keep
    counter <- sum(sc_object_cut$neuron_level==cell)
    #print(cell)
    #print(counter)
    if(counter >= 30){
      keep_list <- sample(x = sc_object_cut$Barcode[sc_object_cut$neuron_level==cell], size = 30, replace = F)
      return(keep_list)
    }
    else { keep_list <- sc_object_cut$Barcode[sc_object_cut$neuron_level==cell]
    return(keep_list)} }) |> unlist()
  
  res <- sc_object_cut[,sc_object_cut$Barcode %in% to_keep]
  
  #print(1)
  
  Idents(object = res) <- 'neuron_level'
  CeNGEN_max_arithMean <- sapply(unique(res$neuron_level), function(cell){
    #print(cell)
    #print(dim(res@assays$RNA@counts[,res$neuron_level==cell]))
    temp_mat <- res@assays$RNA@layers$counts[,res$neuron_level==cell]
    rownames(temp_mat) <- rownames(res@assays$RNA@features@.Data)
    temp_mat <- sweep(temp_mat, 2, colSums(temp_mat), '/')
    return(Matrix::rowMeans(temp_mat))
  })
  CeNGEN_max_arithMean <- CeNGEN_max_arithMean * 1000000
  
  
  CeNGEN_max_arithMean <- CeNGEN_max_arithMean[Matrix::rowSums(CeNGEN_max_arithMean > 10) > 0,]
  
  CeNGEN_max_arithMean_cv <- apply(CeNGEN_max_arithMean, 1, sd)
  CeNGEN_max_arithMean_cv <- CeNGEN_max_arithMean_cv/Matrix::rowMeans(CeNGEN_max_arithMean)
  
  
  CeNGEN_max_arithMean_cut <- CeNGEN_max_arithMean[CeNGEN_max_arithMean_cv > 3,]
  
  #sum(CeNGEN_12max_arithMean_cv > 10, na.rm = T)
  CeNGEN_max_arithMean_CPM.df <- data.frame(CeNGEN_max_arithMean_cut)
  CeNGEN_max_arithMean_CPM.df <- CeNGEN_max_arithMean_CPM.df[order(rownames(CeNGEN_max_arithMean_CPM.df)), 
                                                             order(colnames(CeNGEN_max_arithMean_CPM.df))]
  
  CeNGEN_max_arithMean_CPM.df$DD <- CeNGEN_max_arithMean_CPM.df$VD_DD
  CeNGEN_max_arithMean_CPM.df$VD <- CeNGEN_max_arithMean_CPM.df$VD_DD
  
  common.genes <- intersect(rownames(CeNGEN_max_arithMean_CPM.df), rownames(bulk_data))
  
  est <- sapply(colnames(bulk_data), function(sample1){
    cell <- str_split_fixed(sample1, 'r', 2)[,1]
    cell_types <- c(cell, non_neuronal_list)
    
    CeNGEN_CPMcut <- CeNGEN_max_arithMean_CPM.df[common.genes, cell_types]
    
    
    
    
    bulk <- as.numeric(bulk_data[common.genes, sample1])
    names(bulk) <- common.genes
    
    cell_types
    
    nnls(b=log1p(as.matrix(bulk)), A=log1p(as.matrix(CeNGEN_CPMcut)))$x
    
  })
  
  rownames(est) <- c('Neuron', non_neuronal_list)
  est <- sweep(est, 2, colSums(est), '/')
  
  return(est)
  
})
NNLS_30_list_log1p[[1]]

NNLS_30_list_log1p[[1]] |> colSums()

NNLS_reduce_log1p <- Reduce('+', NNLS_30_list_log1p)/length(NNLS_30_list_log1p)
NNLS_reduce_log1p
write.table(NNLS_reduce_log1p, 'Data_out/NNLS_average_across_100_bootstraps.log1p.30Cells.231211.tsv', sep = '\t')

NNLS_reduce_log1p <- read.table('Data_out/NNLS_average_across_100_bootstraps.log1p.30Cells.231211.tsv')

colSums(NNLS_reduce_log1p)

NNLS_reduce_log1p[,1:5]

sc_object_cut$Cell.type


sc_size <- sapply(unique(sc_object_cut$Cell.type), function(cell){
  return(sum(sc_object_cut$Cell.type==cell))
})
sc_size['DD'] <- sc_size['VD_DD']
sc_size['VD'] <- sc_size['VD_DD']

sc_size <- sc_size[order(names(sc_size))]
sc_size <- sc_size[-length(sc_size)]

proportions[,1:5]


data.frame(sc_size = log10(sc_size[str_split_fixed(colnames(NNLS_reduce_log1p), 'r', 2)[,1]]),
           Neuron = unlist(NNLS_reduce_log1p['Neuron',])) %>%
  ggplot() + 
  theme_classic(base_size = 20) + 
  geom_point(aes(x = sc_size, y = Neuron), color = 'black', size = 4) +
  geom_smooth(aes(x = sc_size, y = Neuron), method = 'lm', se = T) +
  xlab('Single Cell Cluster Size, Log10') +  ylim(0,1)

downsampled_lm <- lm(Neuron ~ sc_size, 
                     data = data.frame(sc_size = log10(sc_size[str_split_fixed(colnames(NNLS_reduce_log1p), 'r', 2)[,1]]),
                                       Neuron = unlist(NNLS_reduce_log1p['Neuron',])))
summary(downsampled_lm)


Heatmap(NNLS_reduce_log1p, cluster_columns = F)

data.frame(sc_size = log10(sc_size[str_split_fixed(colnames(proportions), 'r', 2)[,1]]),
           Neuron = unlist(proportions['Neuron',])) %>%
  ggplot() + 
  theme_classic(base_size = 20) + 
  geom_point(aes(x = sc_size, y = Neuron), color = 'black', size = 4) +
  geom_smooth(aes(x = sc_size, y = Neuron), method = 'lm', se = F) +
  xlab('Single Cell Cluster Size, Log10') +  ylim(0,1)
