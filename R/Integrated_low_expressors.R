###~~ protein_coding_genes detected somewhere, but not in certain neurons, raw bulk counts ----


## packages ----
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pbapply)
library(tidyverse)
library(wbData)
library(cengenDataSC)
library(nlstools)

### load data, normalized ----

integrated_data <- read.csv('Data_out/Integrated_thresholded/Integrated_bsn12_cpm_unthresholded_070625.csv')


### load single cell ncell reference

ncells <- read.csv('references/sc_size_032322.csv')

### load Alexis single cell thresholds
ws289 <- data.frame(wbData::wb_load_gene_ids(289))
ws281 <- data.frame(wbData::wb_load_gene_ids(281))
rownames(ws289) <- ws289$gene_id
rownames(ws281) <- ws281$gene_id
ws289_prc <- ws289[ws289$biotype=='protein_coding_gene',]
ws281_prc <- ws281[ws281$biotype=='protein_coding_gene',]
cengen_sc_2_bulk <- cengenDataSC::cengen_sc_2_bulk
sc_gene_thresholds <- readRDS('references/211028_genes_categorized_by_pattern.rds')

sum(sapply(sc_gene_thresholds, length))
unlist(sc_gene_thresholds) |> unique() |> length()
nondetected_threshold_genes <- sc_gene_thresholds[['nondetected']]
nondetected_threshold_genes <- intersect(nondetected_threshold_genes, rownames(ws289_prc))
length(nondetected_threshold_genes)

CeNGEN_TPM <- cengenDataSC::cengen_TPM_bulk



not_detected_Seurat <- setdiff(rownames(integrated_data), rownames(CeNGEN_TPM))
not_detected_Seurat <- intersect(not_detected_Seurat, rownames(ws289_prc))
length(not_detected_Seurat)

length(setdiff(not_detected_Seurat, nondetected_threshold_genes))
length(nondetected_threshold_genes)

sum(nondetected_threshold_genes %in% rownames(cengen_sc_2_bulk))

sum(colSums(cengen_sc_2_bulk) == 0)

table(c(not_detected_Seurat, nondetected_threshold_genes)) |> names() |> length()

cengen_sc_2_bulk['WBGene00005540',]


### subset to the genes of interest given either criteria and combine the matrices ----

all_nondetect <- unique(c(not_detected_Seurat, nondetected_threshold_genes)) |> sort()

all_nondetect <- all_nondetect[!all_nondetect %in% rownames(CeNGEN_TPM)]

sum(all_nondetect %in% rownames(CeNGEN_TPM))/length(all_nondetect)
length(all_nondetect)

all_nondetect |> clipr::write_clip()

aggr_nondetected_all <- integrated_data[all_nondetect,]


####


integrated_data_subset_threshold <- (aggr_nondetected_all > 0.13325) * 1
integrated_data_subset_threshold <- data.frame(integrated_data_subset_threshold)
colSums(integrated_data_subset_threshold)

new_genes <- list()

for (col_name in colnames(integrated_data_subset_threshold)) {
  # Get rownames where the column value is 1
  rows_with_1 <- rownames(integrated_data_subset_threshold)[integrated_data_subset_threshold[[col_name]] == 1]
  
  # Create summary for this column
  new_genes[[col_name]] <- rows_with_1
}



### new genes per cell plot


GT_thresholded_FDR_M.df <- data.frame(protein_genes = sapply(new_genes, length),
                                      cells = names(new_genes))

GT_thresholded_FDR_M.df %>% ggplot(data = ., aes(x = cells, y = protein_genes)) + geom_col() +
  #scale_y_continuous(transform = 'log10') +
  xlab('Cell Type') + ylab('Bulk Genes undectected in specific single cell clusters') +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text = element_text(face = 'bold', color = 'black'))




### single cell cluster size


sc_size <- ncells$x
names(sc_size) <- rownames(ncells)



### modeling decay of number of new genes detected by cluster size ----

new_genes.vs.size.df <- data.frame(row.names = names(new_genes),
                                   raw_protein_genes = sapply(new_genes,length),
                                   sc_size = sc_size[names(new_genes)])

new_genes.vs.size.df |>
  ggplot() +
  geom_point(aes(x = log10(sc_size), raw_protein_genes), size = 5) +
  geom_smooth(aes(x = log10(sc_size), raw_protein_genes), method = 'lm') 


model <- nls(raw_protein_genes ~ a * sc_size^b, 
             data = new_genes.vs.size.df,
             start = list(a = 100, b = -0.5))

summary(model)
nls_new_genes_decay <- nls(raw_protein_genes ~ m + (M-m)*exp(-sc_size/alpha),
                           data = new_genes.vs.size.df, algorithm = 'port',
                          
                           start = list(M = 2, alpha = 5, m = 1))

x_range <- seq(32, 
               max(new_genes.vs.size.df$sc_size), 
               length.out = 52)
pred_df <- data.frame(sc_size = x_range)
pred_df$fit <- predict(model, newdata = pred_df)


boots_nls_new_genes_decay <- nlsBoot(model, niter=2000)
boots_nls_new_genes_decay_conf <- 
  nlsBootPredict(boots_nls_new_genes_decay, interval = 'confidence', newdata = pred_df) |> data.frame()
colnames(boots_nls_new_genes_decay_conf) <- c('Median', 'percentile_025', 'percentile_975')

# new_genes.df <- cbind(new_genes.vs.size.df[order(new_genes.vs.size.df$sc_size),], boots_nls_new_genes_decay_conf)

boots_nls_new_genes_decay_conf$sc_size <- pred_df$sc_size
new_genes.vs.size.df$cell = rownames(new_genes.vs.size.df)
ggplot() +
  theme_classic(base_size = 25) +
  geom_point(data = new_genes.vs.size.df,
             aes(x=sc_size, y = raw_protein_genes), size = 3) +
  geom_line(data = boots_nls_new_genes_decay_conf,
            aes(x=sc_size, y= Median), color = "red3", linetype = "dashed") +
  geom_line(data = boots_nls_new_genes_decay_conf,
            aes(x=sc_size, y= percentile_025), color = "blue3", linetype = "dashed") +
  geom_line(data = boots_nls_new_genes_decay_conf,
            aes(x=sc_size, y= percentile_975), color = "blue3", linetype = "dashed") +
  ggrepel::geom_text_repel(data = new_genes.vs.size.df,
                           aes(x=sc_size, y = raw_protein_genes, label = cell), max.overlaps = 2) +
  xlab("Single Cell Cluster size") +
  ylab("'new' protein coding genes found in bulk")



