

### libraries

library(ggplot2)
library(ggrepel)
library(dplyr)
library(wbData)
library(edgeR)
library(pbapply)
library(reshape)
library(stringr)
library(wbData)


#############



neuron_metadata <- read.csv('references/Hammarlund_Kaitlyn_Alec_neuron_annotation_070921.csv')

colnames(neuron_metadata) <- gsub('\\.\\.', '_', colnames(neuron_metadata))
colnames(neuron_metadata) <- gsub('\\.', '_', colnames(neuron_metadata))


neuron_metadata$Modality <- sapply(neuron_metadata$Neuron, function(x){
  
  if(neuron_metadata[neuron_metadata$Neuron==x, "Modality_Motor"] == 1){
    return('Motor')
  } else{ if(neuron_metadata[neuron_metadata$Neuron==x, "Modality_Sensory"] == 1){return('Sensory')}else{ if(neuron_metadata[neuron_metadata$Neuron==x, "Modality_Interneuron"] == 1){return('Interneuron')} else{return('Unknown')}}}
  
})

neuron_metadata[neuron_metadata$Neuron == 'IL1', "Modality"] <- 'Sensory'

#### bulk data read in ----


bulk_data <- read.table('Data/Bulk_data_bsn12_231211.tsv', header = T, sep = '\t')

bulk_data[1:10,1:10]

colnames(bulk_data)

## cut reference samples and outliers

bulk_data$RICr133 <- NULL
bulk_data$PVMr122 <- NULL




### remove rRNA, miRNA, and piRNA genes

ref_genes <- wbData::wb_load_gene_ids('289')

ref_genes$biotype |> unique() |> sort()

remove_list <- c('antisense_lncRNA_gene', 'gene', 'miRNA_gene', 'rRNA_gene', 'piRNA_gene')

biotypes_keep <- ref_genes |>
  dplyr::filter(!biotype %in% remove_list) |>
  pull(gene_id)


## which cell types have > 1 replicate?

samples <- colnames(bulk_data)
samples_cell_types <- do.call(rbind, strsplit(colnames(bulk_data), 'r', 2))[,1]
cell_types <- unique(samples_cell_types)



### running edgeR

bulk_data <- bulk_data[rownames(bulk_data) %in% biotypes_keep,]

sum(is.na(bulk_data))

genes_keep <- rowSums(bulk_data > 5) > 2

sum(genes_keep)

# bulk_data[rowSums(bulk_data > 5) > 2,] |> dim()

bulk_data <- bulk_data[genes_keep,]


#bulk_reps2keep_cell_types <- str_split_fixed(bulk_reps2keep, 'r', 2)[,1]

bulk_data |> colnames()

raw_edgeR <- edgeR::DGEList(counts = bulk_data, group = samples_cell_types)
raw_edgeR <- calcNormFactors(raw_edgeR, method = 'TMM')

bulk_tmm <- cpm(raw_edgeR, normalized.lib.sizes = T)

aggr_raw_TMM <- bulk_tmm
colnames(aggr_raw_TMM) <-str_split_fixed(colnames(aggr_raw_TMM),"r",2)[,1]
aggr_raw_TMM <- data.frame(vapply(unique(colnames(aggr_raw_TMM)), function(x)
  rowMeans(aggr_raw_TMM[,colnames(aggr_raw_TMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_TMM)) ))
dim(aggr_raw_TMM)


bulk_tmm_log1p <- bulk_tmm |> log1p()
aggr_raw_TMM_log1p <- aggr_raw_TMM |> log1p()

cv <- apply(bulk_tmm_log1p, 1, function(x){sd(x)/mean(x)})

hvg <- cv[order(cv, decreasing = T)][1:10000] |> names()

bulk_pca <- bulk_tmm_log1p[hvg,] |> t() |> prcomp()

bulk_pca_df <- bulk_pca$x[,1:5] |> data.frame()

bulk_pca_df$sample <- rownames(bulk_pca_df)
bulk_pca_df$neuron <- str_split_fixed(bulk_pca_df$sample, 'r', 2)[,1]

bulk_pca_df$modality <- sapply(bulk_pca_df$neuron, function(x){neuron_metadata[neuron_metadata$Neuron==x, "Modality"]})


ggplot(bulk_pca_df) + geom_label(aes(x = PC1, y = PC2, label = sample, fill = neuron), alpha = 0.8) +
  theme_classic(base_size = 20) +
  theme(legend.position = '', axis.title = element_text(colour = 'black', face = 'bold'))
ggsave('bulk RNA-seq pca neuron label 240527.pdf', width = 10, height = 7)

ggplot(bulk_pca_df) + geom_label(aes(x = PC1, y = PC2, label = sample, fill = modality), alpha = 0.8) +
  theme_classic(base_size = 20) +
  theme(legend.position = '')
ggsave('bulk RNA-seq pca modality label 240804.pdf', width = 10, height = 7)


bulk_pca_df_mean <- bulk_pca_df |> group_by(neuron) |>
  summarize(PC1 = mean(PC1), PC2 = mean(PC2),
            neuron = unique(neuron),
            modality = unique(modality))

ggplot(bulk_pca_df_mean) + geom_label(aes(x = PC1, y = PC2, label = neuron, fill = modality), alpha = 0.8) +
  theme_classic(base_size = 20) +
  theme(legend.position = '')


cv_mean <- apply(aggr_raw_TMM_log1p, 1, function(x){sd(x)/mean(x)})

hvg_mean <- cv_mean[order(cv_mean, decreasing = T)][1:20000] |> names()

aggr_bulk_pca <- aggr_raw_TMM_log1p[hvg_mean,] |> t() |> prcomp()

aggr_bulk_pca_df <- aggr_bulk_pca$x[,1:5] |> data.frame()

aggr_bulk_pca_df$neuron <- rownames(aggr_bulk_pca_df)

aggr_bulk_pca_df$modality <- sapply(aggr_bulk_pca_df$neuron, function(x){neuron_metadata[neuron_metadata$Neuron==x, "Modality"]})




ggplot(aggr_bulk_pca_df, aes(x = PC1, y = -PC2, label = neuron, fill = modality, color = modality)) + 
  #geom_point() +
  geom_label_repel(color = 'black', alpha = 0.8) +
  theme_classic(base_size = 20) +
  theme(legend.position = '')
ggsave('figures/bulk RNA-seq pca cell average modality label 241106.pdf', width = 10, height = 7)


