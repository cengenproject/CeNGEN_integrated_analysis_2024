### detection of single cell genes in bulk
library(dendextend)
library(stringr)
library(ComplexHeatmap)
library(pbapply)
library(edgeR)

### load data

strict <- cengenDataSC::cengen_sc_4_bulk
strict[strict>0] = 1


bulk_data <- read.table('Data/Bulk_data_bsn12_231211.tsv.gz')

bulk_meta <- read.table('Data/bulk_bsn12_metadata.tsv.gz', sep = '\t')
bulk_meta <- bulk_meta[rownames(bulk_data),]


bulk_data_pk <- (bulk_data/bulk_meta$Length) * 1000


bulk_raw_GeTMM <- DGEList(bulk_data_pk, group = str_split_fixed(colnames(bulk_data_pk), 'r', 2)[,1])
bulk_raw_GeTMM <- calcNormFactors(bulk_raw_GeTMM)
bulk_raw_GeTMM <- cpm(bulk_raw_GeTMM, normalized.lib.sizes = T)

aggr_raw_GeTMM <- bulk_raw_GeTMM
colnames(aggr_raw_GeTMM) <-str_split_fixed(colnames(aggr_raw_GeTMM),"r",2)[,1]
aggr_raw_GeTMM <- data.frame(vapply(unique(colnames(aggr_raw_GeTMM)), function(x)
  rowMeans(aggr_raw_GeTMM[,colnames(aggr_raw_GeTMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_GeTMM)) ))



genes_expressed_sc_4_list <- sapply(colnames(strict), function(cell){
  rownames(strict[strict[,cell] > 0,])
})
strict_bulk_detection <- sapply(colnames(bulk_raw_GeTMM), function(cell){
  samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]
  common.genes <- intersect(rownames(aggr_raw_GeTMM), genes_expressed_sc_4_list[[cell]])
  
  bulk <- bulk_raw_GeTMM[common.genes, samples_cell_type %in% c(cell)]
  
  bulk_bin <- bulk > 5
  bulk_bin_ave <- rowMeans(bulk_bin)
  return(sum(bulk_bin_ave > 0.65)/length(common.genes))
  
  
})



strict_TPM <- cengenDataSC::cengen_TPM_bulk[rownames(strict),colnames(strict)] * strict


strict_bulk_corr_pairwise <- pbsapply(colnames(aggr_raw_GeTMM), function(cell){

  
  #common.genes <- intersect(rownames(aggr_raw_GeTMM), rownames(strict_TPM))
  
  all_corr <- sapply(colnames(aggr_raw_GeTMM), function(cell2){
    
    common.genes <- intersect(rownames(aggr_raw_GeTMM), genes_expressed_sc_4_list[[cell2]])
    
    
    Bulk <- aggr_raw_GeTMM[common.genes,cell] |> log1p()
    SC <- strict_TPM[common.genes,cell2] |> log1p()
    
    return(cor(Bulk, SC, method = 'pearson'))
  })
  
  return(all_corr)
  
})


col_fun_ <- circlize::colorRamp2(colors = c('white', '#Cc0202'), breaks = c(min(unlist(strict_bulk_corr_pairwise)),
                                                                       max(unlist(strict_bulk_corr_pairwise))))



col_dend = as.dendrogram(hclust(dist(t(strict_bulk_corr_pairwise))))

pdf('figures/Supplementary figure 1 correlation heatmap.pdf', height = 12, width = 12)
Heatmap(strict_bulk_corr_pairwise, cluster_rows = col_dend, cluster_columns = col_dend, 
        row_names_side = 'left',
        col = col_fun_,
        show_row_dend = F, show_column_dend = F, name = ' ',
        row_title = 'Single Cell', column_title = 'Bulk', column_title_side = 'bottom',
        column_title_gp = gpar(fontsize = 30, fontface = "bold"), 
        row_title_gp = gpar(fontsize = 30, fontface = "bold"),
        column_names_gp = gpar(fontsize = 15, fontface = "bold"), 
        row_names_gp = gpar(fontsize = 15, fontface = "bold"))
dev.off()


