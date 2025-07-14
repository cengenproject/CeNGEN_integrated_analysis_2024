###~~ protein_coding_genes detected somewhere, but not in certain neurons, raw bulk counts ----


# find a list of genes that are expressed in single cell somewhere -> main_list

## per cell type, list genes expressed cell_expressed, then use that to get the setdiff(main_list, cell_expressed) -> non_expressed
### subset it to exclude the non-neuronal ONLY genes
## in the bulk samples, find how many genes are expressed above the overall 5% FPR rate (non_neuronal) and ~14% FDR rate for neuronal ground truth
## 


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

integrated_data <- read.table('Data/bsn12_subtracted_integrated_propadjust_070625.tsv.gz')


### load single cell ncell reference

ncells <- read.csv('references/sc_size_032322.csv')

### load Alexis single cell thresholds
ws289 <- data.frame(wbData::wb_load_gene_ids(289))
rownames(ws289) <- ws289$gene_id
ws289_prc <- ws289[ws289$biotype=='protein_coding_gene',]
protein_coding_genes <- ws289_prc$gene_id



#### 


integrated_data_medium <- integrated_data
integrated_data_medium[integrated_data_medium < 0.23774] <- 0


integrated_expressed_med_list <- sapply(colnames(integrated_data_medium), function(x){
  
  rownames(integrated_data_medium[integrated_data_medium[,x] > 0,])
  
})


cengen_sc_2_bulk <- data.frame(cengen_sc_2_bulk)
cengen_sc_2_bulk$VD_DD <- cengen_sc_2_bulk$VD


sc_2_expr_list <- sapply(colnames(integrated_data_medium), function(x){

  rownames(cengen_sc_2_bulk[cengen_sc_2_bulk[[x]] > 0,])
  
})


sapply(sc_2_expr_list, length)

integrated_expressed_sums <- sapply(integrated_expressed_med_list, length) |> {\(x) x[x>0]}()




medium_threshold_detected <- sapply(colnames(integrated_data_medium), function(x){
  print(x)
  integrated_sum = sum(integrated_data_medium[,x] > 0)
  
  dynamic_props_sum = sum(cengen_sc_2_bulk[,x] > 0)
  
  return(c('cell' = x, 'integrated_sum' = integrated_sum, 'dynamic_props_sum' = dynamic_props_sum))
  
}) |> t() |> data.frame() |>
  mutate(integrated_sum = as.numeric(integrated_sum),
         dynamic_props_sum = as.numeric(dynamic_props_sum),
         integrated_sc_ratio = (integrated_sum/dynamic_props_sum))




medium_threshold_detected_melt <- medium_threshold_detected |> reshape::melt()

cells_ <- unique(medium_threshold_detected_melt$cell)


medium_threshold_detected_melt$ncell <- ncells[medium_threshold_detected_melt$cell,]
medium_threshold_detected_melt$cell_order <- factor(medium_threshold_detected_melt$cell, 
                                                    levels = cells_[order(ncells[cells_,])])

medium_threshold_detected_melt |> filter(variable != 'integrated_sc_ratio') |>
  ggplot() +
  geom_col(aes(x = cell, y = value, fill = variable), position = 'dodge2') +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(face = 'bold', color = 'black'))
ggsave('figures/Dynamic_props_specific_enrichment/Medium_threshold_detected_genes_070625.pdf',
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
med_disagreement_sc$IL1 |> clipr::write_clip()


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
med_disagreement_int$IL1 |> clipr::write_clip()



med_agreement <- sapply(names(sc_2_expr_list), function(x){
  
  sc <- sc_2_expr_list[[x]]
  inte <- integrated_expressed_med_list[[x]]
  
  return(intersect(inte, sc))
  
})





### writing supplementary tables

max_length <- max(sapply(med_disagreement_int, length))
padded_list <- lapply(med_disagreement_int, function(x) {
  pad_length <- max_length - length(x)
  x <- c(x, rep('', pad_length))
  return(x)
})
df <- data.frame(padded_list)
df |> View()
write.csv(df, 'Data_out/integrated_exclusive_genes.csv')



