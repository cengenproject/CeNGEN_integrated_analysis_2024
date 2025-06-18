### generate integrated agggregate profiles for each neuron
library(dplyr)
library(edgeR)
library(stringr)

### function for integrating ----

integrate_geometricMean_biorep <- function(bulk_replicates, 
                                           single_cell_replicates, 
                                           pseudocount, 
                                           bulk_sep = 'r', 
                                           single_cell_sep = '__'){
  
  common.genes <- intersect(rownames(bulk_replicates), rownames(single_cell_replicates))
  
  bulk_replicates <- bulk_replicates[common.genes,]
  single_cell_replicates <- single_cell_replicates[common.genes,]
  
  bulk_cell_types <- str_split_fixed(colnames(bulk_replicates), bulk_sep, 2)[,1]
  cells_bio_rep <- str_split_fixed(colnames(single_cell_replicates), single_cell_sep, 2)[,1]
  
  index_list_sc <- numeric()
  index_list_bulk <- numeric()
  for(g in unique(bulk_cell_types)){
    rep_totals <- c(sum(bulk_cell_types==g), sum(cells_bio_rep==g))
    rep_min <- min(rep_totals)
    #print(g)
    #print(sum(bulk_cell_types==g))
    #print(rep_min)
    indices_sc <- sample(which(cells_bio_rep==g), size = rep_min)
    index_list_sc <- c(index_list_sc, indices_sc)
    indices_bulk <- sample(which(bulk_cell_types==g), size = rep_min)
    index_list_bulk <- c(index_list_bulk, indices_bulk)
    
  }
  #print(index_list_sc)
  #print(index_list_bulk)
  
  #print('set indexes')
  bulk_replicates_match <- bulk_replicates[,index_list_bulk]
  
  single_cell_replicates_match <- single_cell_replicates[,index_list_sc]
  single_cell_replicates_match <- sweep(single_cell_replicates_match, 2, colSums(single_cell_replicates_match, na.rm = T), '/') * colSums(bulk_replicates_match, na.rm = T)
  
  #bulk_replicates_match <- bulk_replicates_match[,order(colnames(bulk_replicates_match))]
  #single_cell_replicates_match <- single_cell_replicates_match[,order(colnames(single_cell_replicates_match))]
  #print(colnames(bulk_replicates_match))
  #print(colnames(single_cell_replicates_match))
  
  sample_level_integration <- exp( ( log(bulk_replicates_match + pseudocount) +
                                       log(single_cell_replicates_match + pseudocount) ) /2 ) - pseudocount
  
  return(sample_level_integration)}



### load in bulk data ----



bulk_subtracted_TMM <- read.table('Data/bsn12_bulk_subtracted_TMM_051624.tsv.gz', sep = '\t')


aggr_subtracted_TMM <- bulk_subtracted_TMM
colnames(aggr_subtracted_TMM) <- str_split_fixed(colnames(aggr_subtracted_TMM),"r",2)[,1]
### combine VD & DD datasets to match single cell annotation
colnames(aggr_subtracted_TMM)[colnames(aggr_subtracted_TMM) %in% c('VD', 'DD')] <- 'VD_DD' 
aggr_subtracted_TMM <- aggr_subtracted_TMM[,order(colnames(aggr_subtracted_TMM))]
aggr_subtracted_TMM <- data.frame(vapply(unique(colnames(aggr_subtracted_TMM)), function(x)
  rowMeans(aggr_subtracted_TMM[,colnames(aggr_subtracted_TMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_subtracted_TMM)) ))


#~~ Single cell ----

# note importing from sc project
prop_by_type <- read.table('Data/SingleCell_proportions_Bulk_annotations.tsv.gz')

colnames(prop_by_type)
## remove non-neuronal profiles
prop_by_type <- prop_by_type[,setdiff(colnames(prop_by_type), c('Intestine',
                                                                'Glia',
                                                                'Pharynx',
                                                                'Excretory',
                                                                'Hypodermis',
                                                                'Rectal_cells',
                                                                'Reproductive',
                                                                'Muscle_mesoderm'))]


#~ Definitions ----

# hard thresholds fixed
low_hard <- 0.02
high_hard <- 0.01


prop_by_type_adjusted <- prop_by_type

prop_by_type_adjusted[apply(prop_by_type_adjusted, 1, min) > high_hard, ] <- 1
prop_by_type_adjusted[apply(prop_by_type_adjusted, 1, max) < low_hard, ] <- 0
mid_range <- (apply(prop_by_type_adjusted, 1, max) > low_hard) & (apply(prop_by_type_adjusted, 1, min) < high_hard)
prop_by_type_adjusted[mid_range, ] <-
  prop_by_type_adjusted[mid_range, ]/apply(prop_by_type_adjusted[mid_range,], 1, max)



sc_biorep <- readRDS('Data/CeNGEN_biorep_aggregate_counts_data_240709.RDS')


celltype_ncell <- sapply(colnames(prop_by_type_adjusted), function(x){
  
  ncell <- sc_biorep$nCells
  
  ncell_x <- ncell[grepl(x, names(ncell))]
  sum(ncell_x)
  
})


sc_counts_from_prop_adjusted <- exp((log(prop_by_type_adjusted/(1-prop_by_type_adjusted))) + log(mean(celltype_ncell)))
sc_counts_from_prop_adjusted[sc_counts_from_prop_adjusted==Inf] <- max(sc_counts_from_prop_adjusted[sc_counts_from_prop_adjusted!=Inf])


common.genes <- intersect(rownames(aggr_subtracted_TMM), rownames(sc_counts_from_prop_adjusted))

neurons <- intersect(colnames(aggr_subtracted_TMM), colnames(sc_counts_from_prop_adjusted))


integrated_data <- apply(abind::abind(sc_counts_from_prop_adjusted[common.genes, neurons],
                                      aggr_subtracted_TMM[common.genes, neurons],
                                      along = 3), c(1,2), function(p){
                                        exp(mean(log(1+p)))-1
                                      })

integrated_data_CPM <- sweep(integrated_data, 
                             MARGIN = 2, 
                             STATS = colSums(integrated_data),
                             FUN = '/') * 1000000

write.table(integrated_data_CPM, 'Data/bsn12_subtracted_integrated_propadjust_071724.tsv',
            sep = '\t', quote = F)
