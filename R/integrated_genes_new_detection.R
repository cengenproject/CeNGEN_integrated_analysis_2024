###~~ protein_coding_genes detected somewhere, but not in certain neurons, raw bulk counts ----


# find a list of genes that are expressed in single cell somewhere -> main_list

## per cell type, list genes expressed cell_expressed, then use that to get the setdiff(main_list, cell_expressed) -> non_expressed
### subset it to exclude the non-neuronal ONLY genes
## in the bulk samples, find how many genes are expressed above the overall 5% FPR rate (non_neuronal) and ~14% FDR rate for neuronal ground truth
## 


## packages ----
library(edgeR)
library(pbapply)
library(tidyverse)
library(wbData)
library(cengenDataSC)
library(nlstools)

### load data, normalized ----

integrated_data <- read.table('Data/bsn12_subtracted_integrated_propadjust_071724.tsv')


### load single cell ncell reference

ncells <- read.csv('references/sc_size_032322.csv')

### load Alexis single cell thresholds
ws289 <- data.frame(wbData::wb_load_gene_ids(289))
rownames(ws289) <- ws289$gene_id
ws289_prc <- ws289[ws289$biotype=='protein_coding_gene',]

sc_gene_thresholds <- readRDS('references/211028_genes_categorized_by_pattern.rds')

nondetected_threshold_genes <- sc_gene_thresholds[['nondetected']]
nondetected_threshold_genes <- intersect(nondetected_threshold_genes, rownames(ws289_prc))


CeNGEN_TPM <- cengenDataSC::cengen_TPM_bulk


not_detected_Seurat <- setdiff(rownames(integrated_data), rownames(CeNGEN_TPM))
not_detected_Seurat <- intersect(not_detected_Seurat, rownames(ws289_prc))
length(not_detected_Seurat)


### subset to the genes of interest given either criteria and combine the matrices ----
aggr_not_detected_Seurat <- integrated_data[not_detected_Seurat,]
aggr_nondetected_threshold_genes <- integrated_data[nondetected_threshold_genes,]
dim(aggr_not_detected_Seurat)
dim(aggr_nondetected_threshold_genes)

aggr_nondetected_all <- rbind(aggr_nondetected_threshold_genes, aggr_not_detected_Seurat)
aggr_nondetected_all <- na.omit(aggr_nondetected_all)
dim(aggr_nondetected_all)

names(sc_gene_thresholds)
length(sc_gene_thresholds[[1]])

protein_coding_genes <- ws289_prc$gene_id

ubiquitous_genes <- sc_gene_thresholds[['ubiquitous']]
nondetected_threshold_genes <- sc_gene_thresholds[['nondetected']]
nondetected_threshold_genes <- intersect(nondetected_threshold_genes, ws289_prc$gene_id)
non_neuronal_genes <- sc_gene_thresholds[['nonneuronal']]

length(nondetected_threshold_genes)

cengen_sc_1_bulk <- cengenDataSC::cengen_sc_1_bulk

hist(log10(rowMeans(CeNGEN_TPM[ubiquitous_genes,])))


cengen_sc_1_bulk <- cengen_sc_1_bulk |> data.frame()

cengen_sc_1_bulk$VD_DD <- cengen_sc_1_bulk$VD

### For each cell type that we have bulk data for: ----
###     make a list of genes that are: 1) Not detected in that cell in single cell at the liberal threshold
###                                    2) Are not labeled "undetected" in all single cell tissues
###                                    3) Are not detected exclusively in non-neuronal single cell clusters
sc_unexpr_1_list <- lapply(colnames(integrated_data), function(cell){
  expr <- rownames(cengen_sc_1_bulk[cengen_sc_1_bulk[,cell] > 0,])
  unexpr <- setdiff(rownames(integrated_data), expr)
  unexpr_protein <- intersect(unexpr, protein_coding_genes)
  unexpr_but_detected <- setdiff(unexpr_protein, nondetected_threshold_genes)
  unexpr_but_not_non_neuronal <- setdiff(unexpr_but_detected, non_neuronal_genes)
  
  return(unexpr_but_not_non_neuronal)
})
names(sc_unexpr_1_list) <- colnames(integrated_data)
sapply(sc_unexpr_1_list, length)


### get a list of the potential genes, no duplicates
all_potential_unexpr <- unique(unlist(sc_unexpr_1_list))
length(all_potential_unexpr)




## calling the genes ----




GT_thresholded_FDR_M <- sapply(colnames(integrated_data), function(cell){
  x1 <- rownames(integrated_data[sc_unexpr_1_list[[cell]],][integrated_data[sc_unexpr_1_list[[cell]], cell] >= 0.23774,])
  length(x1)
})

GT_thresholded_FDR_M.gene.list <- sapply(colnames(integrated_data), function(cell){
  rownames(integrated_data[sc_unexpr_1_list[[cell]],][integrated_data[sc_unexpr_1_list[[cell]], cell] >= 0.23774,])
  
})


integrated_data_medium <- integrated_data
integrated_data_medium[integrated_data_medium < 0.23774] <- 0


integrated_expressed_med_list <- sapply(colnames(integrated_data_medium), function(x){
  
  rownames(integrated_data_medium[integrated_data_medium[,x] > 0,])
  
})


sc_1_expr_list <- sapply(colnames(integrated_data_medium), function(x){
  
  rownames(cengen_sc_1_bulk[cengen_sc_1_bulk[[x]] > 0,])
  
})
cengen_sc_2_bulk <- cengen_sc_2_bulk |> data.frame()
cengen_sc_2_bulk$VD_DD <- cengen_sc_2_bulk$VD

sc_2_expr_list <- sapply(colnames(integrated_data_medium), function(x){
  
  rownames(cengen_sc_2_bulk[cengen_sc_2_bulk[[x]] > 0,])
  
})



rowSums(cengen_sc_1_bulk)

sapply(sc_1_expr_list, length)
sapply(sc_2_expr_list, length)

integrated_expressed_sums <- sapply(integrated_expressed_med_list, length) |> {\(x) x[x>0]}()




medium_threshold_detected <- sapply(colnames(integrated_data_medium), function(x){
  
  integrated_sum = sum(integrated_data_medium[,x] > 0)
  
  dynamic_props_sum = sum(cengen_sc_2_bulk[,x] > 0)
  
  return(c('cell' = x, 'integrated_sum' = integrated_sum, 'dynamic_props_sum' = dynamic_props_sum))
  
}) |> t() |> data.frame() |>
  mutate(integrated_sum = as.numeric(integrated_sum),
         dynamic_props_sum = as.numeric(dynamic_props_sum),
         integrated_sc_ratio = (integrated_sum/dynamic_props_sum))




medium_threshold_detected_melt <- medium_threshold_detected |> reshape::melt()

cells_ <- unique(medium_threshold_detected_melt$cell)

cells_[order(ncells[cells_])]


medium_threshold_detected_melt$ncell <- ncells[medium_threshold_detected_melt$cell,]
medium_threshold_detected_melt$cell_order <- factor(medium_threshold_detected_melt$cell, 
                                                    levels = cells_[order(ncells[cells_,])])
medium_threshold_detected_melt |> filter(variable == 'integrated_sc_ratio') |>
  ggplot() +
  geom_col(aes(x = cell_order, y = log2(value), fill = variable), position = 'dodge2')

medium_threshold_detected_melt |> filter(variable != 'integrated_sc_ratio') |>
  ggplot() +
  geom_col(aes(x = cell_order, y = value, fill = variable), position = 'dodge2') +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(face = 'bold', color = 'black'))
ggsave('figures/Dynamic_props_specific_enrichment/Medium_threshold_detected_genes_240803.pdf',
       height = 5, width = 12)




med_disagreement_sc <- sapply(names(sc_2_expr_list), function(x){
  
  sc <- sc_2_expr_list[[x]]
  inte <- integrated_expressed_med_list[[x]]
  
  return(sc[! sc %in% inte])
  
})

names(med_disagreement_sc) <- names(sc_2_expr_list)
sapply(med_disagreement_sc, length)
med_disagreement_sc$ADL |> clipr::write_clip()
med_disagreement_sc$DVC |> clipr::write_clip()
med_disagreement_sc$VB |> clipr::write_clip()
med_disagreement_sc$AWB |> clipr::write_clip()


med_disagreement_int <- sapply(names(sc_2_expr_list), function(x){
  
  sc <- sc_2_expr_list[[x]]
  inte <- integrated_expressed_med_list[[x]]
  
  return(inte[! inte %in% sc])
  
})
sapply(med_disagreement_int, length) |> sort()

med_disagreement_int$ADL |> clipr::write_clip()
med_disagreement_int$DVC |> clipr::write_clip()
med_disagreement_int$VB |> clipr::write_clip()
med_disagreement_int$AWB |> clipr::write_clip()


plot(sapply(sc_2_expr_list, length)[names(integrated_expressed_sums)], integrated_expressed_sums)
abline(a = 0, b = 1)


names(GT_thresholded_FDR_M.gene.list) <- colnames(integrated_data)


apply(integrated_data, 2, function(x){sum(x>0.23774)})


GT_thresholded_FDR_M.df <- data.frame(protein_genes = as.numeric(GT_thresholded_FDR_M),
                                      cells = names(GT_thresholded_FDR_M))

GT_thresholded_FDR_M.gene.list$OLQ %>% clipr::write_clip(.)

sapply(GT_thresholded_FDR_M.gene.list, length)


cell_names <- rownames(ncells)
sc_size <- ncells$x
names(sc_size) <- cell_names

sc_size[names(GT_thresholded_FDR_M.gene.list)][order(sc_size[names(GT_thresholded_FDR_M.gene.list)])]

IL1_show_TMM <- integrated_data[GT_thresholded_FDR_M.gene.list$IL1 ,]


rownames(IL1_show_TMM) <- wbData::i2s(GT_thresholded_FDR_M.gene.list$IL1 , wb_load_gene_ids('277'))
#View(IL1_show_TMM[,c('AFD', 'I5', 'IL1', 'IL2')])


conf.model <- lm(GT_thresholded_FDR_M ~ 1)
confint(conf.model)
mean(GT_thresholded_FDR_M) - confint(conf.model)[[1]]

##### bar plot of "new" genes ----
GT_thresholded_FDR_M.df %>% ggplot(data = ., aes(x = cells, y = protein_genes)) + geom_col() +
  #scale_y_continuous(transform = 'log10') +
  xlab('Cell Type') + ylab('Bulk Genes undectected in specific single cell clusters') +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text = element_text(face = 'bold', color = 'black'))

ggsave('figures/ncRNA/supp-low expressed genes/New_LowExpressed_Protein_coding_genes_240803.pdf', width = 9, height = 6)

GT_thresholded_FDR_M.df$protein_genes |> quantile()

mean(GT_thresholded_FDR_M.df$protein_genes)
min(GT_thresholded_FDR_M.df$protein_genes)
max(GT_thresholded_FDR_M.df$protein_genes)

sapply(1:1000, function(x){
  set.seed(x)
  t_ <- GT_thresholded_FDR_M.df |> dplyr::sample_n( size = nrow(GT_thresholded_FDR_M.df), replace = T)
  mean(t_$protein_genes)
  
}) |> quantile(c(0.025, 0.975))


new_genes.vs.size.df <- data.frame(row.names = GT_thresholded_FDR_M.df$cells,
                                   raw_protein_genes = GT_thresholded_FDR_M.df$protein_genes,
                                   sc_size = sc_size[GT_thresholded_FDR_M.df$cells])
new_genes.vs.size.df['OLL',]

new_genes.vs.size.df %>%
  ggplot() +
  geom_point(aes(x = log10(sc_size), log10(raw_protein_genes)), size = 5) +
  geom_smooth(aes(x = log10(sc_size), log10(raw_protein_genes)), method = 'lm') 


library(nlstools)
nls_new_genes_decay <- nls(raw_protein_genes ~ m + (M-m)*exp(-sc_size/alpha), data = new_genes.vs.size.df,
                           list(M = 200, alpha = 200, m = 9))


boots_nls_new_genes_decay <- nlsBoot(nls_new_genes_decay, niter=2000)
boots_nls_new_genes_decay_conf <- nlsBootPredict(boots_nls_new_genes_decay, interval = 'confidence')


new_genes.df <- cbind(new_genes.vs.size.df, boots_nls_new_genes_decay_conf)

new_genes.df[order(new_genes.df$sc_size),]
new_genes.df$cell <- rownames(new_genes.df)
ggplot(new_genes.df) +
  theme_classic(base_size = 25) +
  geom_point(aes(x=sc_size, y = raw_protein_genes), size = 3) +
  geom_line(aes(x=sc_size, y= Median), color = "red3", linetype = "dashed") +
  geom_line(aes(x=sc_size, y= `2.5%`), color = "blue3", linetype = "dashed") +
  geom_line(aes(x=sc_size, y= `97.5%`), color = "blue3", linetype = "dashed") +
  ggrepel::geom_text_repel(aes(x=sc_size, y = raw_protein_genes, label = cell), max.overlaps = 2) +
  xlab("Single Cell Cluster size") +
  ylab("'new' protein coding genes found in bulk")
ggsave('figures/ncRNA/supp-low expressed genes/new_genes_versus_sc_cluster_size_240803.pdf')

new_genes_decay




