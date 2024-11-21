

library(Prop2CountR)
library(Seurat)


### load single cell counts data

sc_object <- readRDS('032224_L4_all_cells_Seurat5.rds')


sc_object <- sc_object[,sc_object@meta.data$Tissue == 'Neuron']
sc_object@meta.data$Cell.type_new <- sc_object@meta.data$Cell.type

sc_object@meta.data$Cell.type_new[grepl('^AWC', sc_object@meta.data$Cell.type_new)] <- 'AWC'
sc_object@meta.data$Cell.type_new[grepl('^DA', sc_object@meta.data$Cell.type_new)] <- 'DA'
sc_object@meta.data$Cell.type_new[grepl('^DB', sc_object@meta.data$Cell.type_new)] <- 'DB'
sc_object@meta.data$Cell.type_new[grepl('^IL2', sc_object@meta.data$Cell.type_new)] <- 'IL2'
sc_object@meta.data$Cell.type_new[grepl('^RMD', sc_object@meta.data$Cell.type_new)] <- 'RMD'
sc_object@meta.data$Cell.type_new[grepl('^RME', sc_object@meta.data$Cell.type_new)] <- 'RME'
sc_object@meta.data$Cell.type_new[grepl('^VB', sc_object@meta.data$Cell.type_new)] <- 'VB'
sc_object@meta.data$Cell.type_new[grepl('^VC', sc_object@meta.data$Cell.type_new)] <- 'VC'

sc_object@meta.data$Cell.type_new |> table() |> names() |> sort()

sc_object@meta.data$Cell.type

get_replicate_counts_proportions_Seurat_ <- function(Seurat_object, idents, min_cells = 10, sep = "__"){
  if (length(idents > 1)) {
    ident <- paste(idents, sep = sep)
    
    identity_matrix <- Seurat_object@meta.data |> tidyr::unite('new', idents, sep = sep)
  }
  else {
    ident <- idents
  } 
  nCells <- table(identity_matrix$new)
  nCells <- nCells[nCells >= min_cells]
  replicates <- sort(names(nCells))
  nCells <- nCells[replicates]
  counts <- pbsapply(replicates, function(x) {
    Matrix::rowSums(Seurat_object@assays$RNA@counts[, identity_matrix$new == 
                                                      x])
  })
  proportions <- pbsapply(replicates, function(x) {
    Matrix::rowSums(Seurat_object@assays$RNA@counts[, identity_matrix$new == 
                                                      x] > 0)/nCells[x]
  })
  return(list(counts = counts, proportions = proportions, nCells = nCells))
}

sc_biorep <- get_replicate_counts_proportions_Seurat_(sc_object, idents = c('Cell.type_new', 'Experiment'), min_cells = 10)


neuron_prop2count <- prop2count(proportions_matrix = sc_biorep$proportions |> 
                                  prevent_infinite(nCell_vector = sc_biorep$nCells),
                                nCell_vector = sc_biorep$nCells)

sc_replicate_number <- stringr::str_split_fixed(colnames(neuron_prop2count), '__', 2)[,1] |> table()
sc_cell_keep <- sc_replicate_number[sc_replicate_number>1]

sc_cell_keep <- sc_cell_keep[!grepl(pattern = 'stressed', x = names(sc_cell_keep))]

sc_replicate_keep <- colnames(neuron_prop2count)[str_split_fixed(colnames(neuron_prop2count), '__', 2)[,1] %in% names(sc_cell_keep)]

neuron_prop2count_use <- neuron_prop2count[,sc_replicate_keep]


neuron_prop2count_use
max(neuron_prop2count_use)


write.table(neuron_prop2count_use, 'Data/sc_biorep_prop2count_050124.tsv', sep = '\t', quote = F)
