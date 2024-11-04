
##devtools::install_github('alecbarrett/LittleBites/')


library(LittleBites)

### ground truth ----

UNN_ground_truth_mtx <- read.table('../references/ubiquitous_non_neuronal/ubituiqtous_and_nonNeuronal_gt_genes_matrix_042222.tsv', sep = '\t')
train_test_split <- readRDS('../references/ubiquitous_non_neuronal/ubituiqtous_and_nonNeuronal_gt_genes_split_042222.rds')

UNN_train <- UNN_ground_truth_mtx[train_test_split$training_genes,]

UNN_train$DD <- UNN_train$VD_DD
UNN_train$VD <- UNN_train$VD_DD

UNN_test <- UNN_ground_truth_mtx[train_test_split$testing_genes,]

UNN_test$DD <- UNN_test$VD_DD
UNN_test$VD <- UNN_test$VD_DD



###

bulk_data <- read.table('Data/Bulk_data_bsn12_231211.tsv', sep = '\t')

sc_object_cut_neuron_aggregates_cpm <- read.table('Data/singleCell_reference.tsv', sep = '\t')

cells <- colnames(sc_object_cut_neuron_aggregates_cpm) |> unique() |> sort()

specificity <- pbapply(sc_object_cut_neuron_aggregates_cpm |> log1p(), 1, LittleBites::Spm)

colnames(sc_object_cut_neuron_aggregates)

contaminants <- c('Excretory', 'Glia', 'Hypodermis', 'Intestine', 'Muscle_mesoderm', 'Pharynx', 'Reproductive')

neurons <- colnames(sc_object_cut_neuron_aggregates_cpm)[!colnames(sc_object_cut_neuron_aggregates_cpm) %in% c(contaminants, 'Rectal_cells')]


cell_types_matrix <- sapply(colnames(bulk_data), function(sample_){
  
  cell <- str_split_fixed(sample_, 'r', 2)[,1]
  
  return(c(cell, contaminants))
  
  
}) |> t()
cell_types_matrix


cell_types_matrix[1:5,]

common.genes <- intersect(rownames(bulk_data),
                          rownames(sc_object_cut_neuron_aggregates_cpm))

sc_object_cut_neuron_aggregates_cpm <- data.frame(sc_object_cut_neuron_aggregates_cpm)

sc_object_cut_neuron_aggregates_cpm$DD <- sc_object_cut_neuron_aggregates_cpm$VD_DD
sc_object_cut_neuron_aggregates_cpm$VD <- sc_object_cut_neuron_aggregates_cpm$VD_DD

bulk_data_use <- bulk_data[common.genes,]
sc_object_cut_neuron_aggregates_cpm_use <- sc_object_cut_neuron_aggregates_cpm[common.genes,]
specificity <- specificity[common.genes]

bulk_subtracted <- subtraction(bulk = bulk_data_use,
                         reference = sc_object_cut_neuron_aggregates_cpm_use,
                         cell_types_matrix = cell_types_matrix, 
                         training_matrix = UNN_train, 
                         specificity_weights = specificity, verbose = F)

write.table(bulk_subtracted, 'Data_out/bsn12_bulk_subtracted_030424.tsv', sep = '\t')



bulk_raw_train_AUROC <- sapply(samples, function(sample_1){
  
  cells_in_use <- cell_types_matrix[sample_1,] |> unlist()
  
  cell <- cells_in_use[1]
  contaminant_tissues <- cells_in_use[2:length(cells_in_use)]
  
  bulk_deconv_target <- bulk_data_use[,sample_1] ## some steps required a dataframe
  names(bulk_deconv_target) <- rownames(bulk)
  
  starting_auc <- calc_bulk_auc_one_sample(bulk_deconv_target,
                                           sample_1,
                                           training_matrix,
                                           training_genes = rownames(training_matrix),
                                           
                                           sep = sep)
  return(starting_auc)
  
  
})

bulk_sub_train_AUROC <- sapply(samples, function(sample_1){
  
  cells_in_use <- cell_types_matrix[sample_1,] |> unlist()
  
  cell <- cells_in_use[1]
  contaminant_tissues <- cells_in_use[2:length(cells_in_use)]
  
  bulk_deconv_target <- bulk_subtracted[,sample_1] ## some steps required a dataframe
  names(bulk_deconv_target) <- rownames(bulk)
  
  starting_auc <- calc_bulk_auc_one_sample(bulk_deconv_target,
                                           sample_1,
                                           training_matrix,
                                           training_genes = rownames(training_matrix),
                                           
                                           sep = sep)
  return(starting_auc)
  
  
})

plot(bulk_raw_train_AUROC, bulk_sub_train_AUROC,
     xlim = c(0.72,1), ylim = c(0.72,1),
     xlab = 'bulk AUROC', ylab = 'subtracted AUROC')
abline(a=0,b=1,col='red')


bulk_raw_test_AUROC <- sapply(samples, function(sample_1){
  
  cells_in_use <- cell_types_matrix[sample_1,] |> unlist()
  
  cell <- cells_in_use[1]
  contaminant_tissues <- cells_in_use[2:length(cells_in_use)]
  
  bulk_deconv_target <- bulk_data_use[,sample_1] ## some steps required a dataframe
  names(bulk_deconv_target) <- rownames(bulk)
  
  starting_auc <- calc_bulk_auc_one_sample(bulk_deconv_target,
                                           sample_1,
                                           UNN_test,
                                           training_genes = rownames(UNN_test),
                                           
                                           sep = sep)
  return(starting_auc)
  
  
})
bulk_sub_test_AUROC <- sapply(samples, function(sample_1){
  
  cells_in_use <- cell_types_matrix[sample_1,] |> unlist()
  
  cell <- cells_in_use[1]
  contaminant_tissues <- cells_in_use[2:length(cells_in_use)]
  
  bulk_deconv_target <- bulk_subtracted[,sample_1] ## some steps required a dataframe
  names(bulk_deconv_target) <- rownames(bulk)
  
  starting_auc <- calc_bulk_auc_one_sample(bulk_deconv_target,
                                           sample_1,
                                           UNN_test,
                                           training_genes = rownames(UNN_test),
                                           
                                           sep = sep)
  return(starting_auc)
  
  
})
plot(bulk_raw_test_AUROC, bulk_sub_test_AUROC,
     xlim = c(0.72,1), ylim = c(0.72,1),
     xlab = 'bulk AUROC', ylab = 'subtracted AUROC')
abline(a=0,b=1,col='red')




bulk_raw_estimate <- sapply(samples, function(sample_1){
  
  cells_in_use <- cell_types_matrix[sample_1,] |> unlist()
  
  cell <- cells_in_use[1]
  contaminant_tissues <- cells_in_use[2:length(cells_in_use)]
  
  bulk_deconv_target <- bulk_data_use[,sample_1] ## some steps required a dataframe
  names(bulk_deconv_target) <- rownames(bulk)
  
  starting_auc <- calc_bulk_auc_one_sample(bulk_deconv_target,
                                           sample_1,
                                           training_matrix,
                                           training_genes = rownames(training_matrix),
                                           
                                           sep = sep)
  reference_tmp <- reference[,cells_in_use]
  reference_tmp <- sweep(reference_tmp, 2, colSums(reference_tmp), '/')
  reference_tmp <- reference_tmp * sum(bulk_deconv_target)
  estimates <- nnls::nnls( A = log1p(as.matrix(reference_tmp * specificity_weights )),
                           b = log1p(bulk_deconv_target * specificity_weights ) )$x
  names(estimates) <- cells_in_use
  estimates <- estimates/sum(estimates)
  
  
  return(estimates)
  
  
})
bulk_sub_estimate <- sapply(samples, function(sample_1){
  
  cells_in_use <- cell_types_matrix[sample_1,] |> unlist()
  
  cell <- cells_in_use[1]
  contaminant_tissues <- cells_in_use[2:length(cells_in_use)]
  
  bulk_deconv_target <- bulk_subtracted[,sample_1] ## some steps required a dataframe
  names(bulk_deconv_target) <- rownames(bulk)
  
  reference_tmp <- reference[,cells_in_use]
  reference_tmp <- sweep(reference_tmp, 2, colSums(reference_tmp), '/')
  reference_tmp <- reference_tmp * sum(bulk_deconv_target)
  estimates <- nnls::nnls( A = log1p(as.matrix(reference_tmp * specificity_weights )),
                           b = log1p(bulk_deconv_target * specificity_weights ) )$x
  names(estimates) <- cells_in_use
  estimates <- estimates/sum(estimates)
  
  
  return(estimates)
  
  
})
plot(bulk_raw_estimate[1,], bulk_sub_estimate[1,],
     xlim = c(0,1), ylim = c(0,1),
     xlab = 'bulk log1p NNLS estimate', ylab = 'subtracted log1p NNLS estimate')
abline(a=0,b=1,col='red')

tmp <- data.frame(orig = bulk_raw_estimate[1,], subtracted = bulk_sub_estimate[1,]) |>
  lm(formula = subtracted~orig, data = _)
tmp |> summary()
lm()
