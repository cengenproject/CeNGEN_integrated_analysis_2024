### libraries

library(dplyr)
library(purrr)
library(ggplot2)
library(bayestestR)
library(pROC)
library(stringr)




### define functions ----
get_tpr <- function(expression, truth, threshold, na.rm = TRUE){
  # True Positive Rate, aka sensitivity, aka recall
  # TPR = TP/(TP+FN) = TP/P
  bin <- expression >= threshold
  return(sum(bin * truth)/sum(truth))
}
get_fpr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Positive Rate
  # FPR = FP/(FP+TN) = FP/N
  bin <- expression >= threshold
  return(sum(bin * (!truth))/sum(!(truth)))
}
get_fdr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Discovery Rate
  # FDR = FP/(FP+TP) = 1 - PPV
  bin <- expression >= threshold
  fdr <- sum(bin * (!truth))/(sum(bin*(!truth)) + sum(bin*truth))
  if(is.nan(fdr))
    fdr <- 0
  return(fdr)
}

#### load data

# save(bmind_neuron_count_TMM, egm_TMM, bulk_subtracted_TMM, bulk_raw_TMM, genes, neurons, testing_matrix, 
#      file = 'AUROC_Check_for_Alec.RData')


#load('../bsn12/AUROC_Check_for_Alec.RData')

testing_matrix <-  read.table('references/ubiquitous_nonNeuronal_testing_matrix.tsv', sep = '\t')




bulk_raw_TMM <- read.table('Data/bsn12_bulk_TMM_051624.tsv.gz', sep = '\t')
bulk_subtracted_TMM <- read.table('Data/bsn12_bulk_subtracted_TMM_051624.tsv.gz', sep = '\t')
bmind_neuron_count_TMM <- read.table('Data/bsn12_bulk_bMIND_TMM_051624.tsv.gz', sep = '\t')
egm_TMM <- read.table('Data/bsn12_bulk_enigma_TMM_051624.tsv.gz', sep = '\t')


aggr_raw_TMM <- bulk_raw_TMM
colnames(aggr_raw_TMM) <-str_split_fixed(colnames(aggr_raw_TMM),"r",2)[,1]
aggr_raw_TMM <- data.frame(vapply(unique(colnames(aggr_raw_TMM)), function(x)
  rowMeans(aggr_raw_TMM[,colnames(aggr_raw_TMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_TMM)) ))


aggr_subtracted_TMM <- bulk_subtracted_TMM
colnames(aggr_subtracted_TMM) <-str_split_fixed(colnames(aggr_subtracted_TMM),"r",2)[,1]
aggr_subtracted_TMM <- data.frame(vapply(unique(colnames(aggr_subtracted_TMM)), function(x)
  rowMeans(aggr_subtracted_TMM[,colnames(aggr_subtracted_TMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_subtracted_TMM)) ))


bmind_ave <- bmind_neuron_count_TMM
colnames(bmind_ave) <-str_split_fixed(colnames(bmind_ave),"r",2)[,1]
bmind_ave <- data.frame(vapply(unique(colnames(bmind_ave)), function(x) 
  rowMeans(bmind_ave[,colnames(bmind_ave)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(bmind_ave)) ))

egm_ave <- egm_TMM
colnames(egm_ave) <-str_split_fixed(colnames(egm_ave),"r",2)[,1]
egm_ave <- data.frame(vapply(unique(colnames(egm_ave)), function(x) 
  rowMeans(egm_ave[,colnames(egm_ave)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(egm_ave)) ))



genes <- rownames(testing_matrix)

genes <- genes[(genes %in% rownames(egm_ave)) & (genes %in% rownames(bmind_ave))]

neurons <- colnames(aggr_subtracted_TMM)

testing_matrix$VD <- testing_matrix$VD_DD
testing_matrix$DD <- testing_matrix$VD_DD

testing_gt <- testing_matrix[genes, neurons]


aggr_raw_TMM_plot <- aggr_raw_TMM[genes, neurons]
aggr_sub_TMM_plot <- aggr_subtracted_TMM[genes, neurons]
egm_ave_plot <- egm_ave[genes, neurons]
bmind_ave_plot <- bmind_ave[genes, neurons]




diags_aggr_raw_ave_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                                  TPR = map_dbl(threshold, ~get_tpr(aggr_raw_TMM_plot, testing_gt, .x)),
                                  FPR = map_dbl(threshold, ~get_fpr(aggr_raw_TMM_plot, testing_gt, .x)),
                                  FDR = map_dbl(threshold, ~get_fdr(aggr_raw_TMM_plot, testing_gt, .x)),
                                  counts = "raw")


diags_aggr_sub_TMM_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                                  TPR = map_dbl(threshold, ~get_tpr(aggr_sub_TMM_plot, testing_gt, .x)),
                                  FPR = map_dbl(threshold, ~get_fpr(aggr_sub_TMM_plot, testing_gt, .x)),
                                  FDR = map_dbl(threshold, ~get_fdr(aggr_sub_TMM_plot, testing_gt, .x)),
                                  counts = "Subtracted Bulk")


diags_aggr_enigma_L2_log_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                                        TPR = map_dbl(threshold, ~get_tpr(egm_ave_plot, testing_gt, .x)),
                                        FPR = map_dbl(threshold, ~get_fpr(egm_ave_plot, testing_gt, .x)),
                                        FDR = map_dbl(threshold, ~get_fdr(egm_ave_plot, testing_gt, .x)),
                                        counts = "ENIGMA")


diags_aggr_bMIND_ave_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                                    TPR = map_dbl(threshold, ~get_tpr(bmind_ave_plot, testing_gt, .x)),
                                    FPR = map_dbl(threshold, ~get_fpr(bmind_ave_plot, testing_gt, .x)),
                                    FDR = map_dbl(threshold, ~get_fdr(bmind_ave_plot, testing_gt, .x)),
                                    counts = "bMIND")

bind_rows(diags_aggr_raw_ave_plot,
          diags_aggr_sub_TMM_plot,
          diags_aggr_enigma_L2_log_plot,
          diags_aggr_bMIND_ave_plot
) |> 
  ggplot(aes(x = FPR, y=TPR, color= counts)) +
  geom_abline(slope = 1) +
  geom_path(linewidth = 2, alpha = 0.7) +
  ggtitle('ROC for ubiquitous and non-neuronal\ntesting genes') +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(color = 'black', face = 'bold'), 
        axis.title = element_text(color = 'black', face = 'bold'),
        title = element_text(color = 'black', face = 'bold'))

ggsave('figures/UN_Testing_ROC_051524.pdf', width = 7, height = 7)

bind_rows(diags_aggr_raw_ave_plot,
          diags_aggr_sub_TMM_plot,
          diags_aggr_enigma_L2_log_plot,
          diags_aggr_bMIND_ave_plot
) |> 
  ggplot(aes(x = 1-FDR, y=TPR, color= counts)) +
  geom_vline(xintercept = .95) +
  geom_path(linewidth = 3, alpha = 0.8) +
  ggtitle('Precision-Recall\nfor ubiquitous and non-neuronal\ntesting genes') +
  theme_classic(base_size = 20) #+ theme(legend.position = '')


raw_roc <- roc(testing_gt |> unlist(),
               aggr_raw_TMM_plot|> log1p() |> unlist())

sub_roc <- roc(testing_gt |> unlist(),
               aggr_sub_TMM_plot |> log1p() |> unlist())

bmind_roc <- roc(testing_gt |> unlist(),
                 bmind_ave_plot |> log1p() |> unlist())

enigma_roc <- roc(testing_gt |> unlist(),
                  egm_ave_plot |> log1p() |> unlist())


roc_list <- list('unaltered bulk' = raw_roc, 'subtracted' = sub_roc, 
              'bMIND' = bmind_roc, 'enigma' = enigma_roc)



sapply(roc_list, function(y){
  x <- ci.auc(y)
  return(c('lower_ci' = x[1], 'mean' = x[2], 'upper_ci' = x[3]))
}) |> t() |> data.frame() |> tibble::rownames_to_column('dataset') |>
  mutate(dataset = dataset |> factor(levels = c('unaltered bulk', 'bMIND', 'enigma', 'subtracted'))) |>
  ggplot() + 
  geom_col(aes(x = dataset, y = mean, fill = dataset), alpha = 0.8) + 
  geom_errorbar(aes(x = dataset, ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  geom_text(aes(x = dataset, y = mean+0.02, label = mean |> round(4))) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black', face = 'bold'), 
    axis.title = element_text(color = 'black', face = 'bold'),
    title = element_text(color = 'black', face = 'bold'))
ggsave('figures/UN_Testing_ROC_barchart_051624.pdf', width = 7, height = 7)




raw_vs_bmind <- roc.test(raw_roc, bmind_roc, method = 'delong')

raw_vs_enigma <- roc.test(raw_roc, enigma_roc, method = 'delong')
raw_vs_subtracted <- roc.test(raw_roc, sub_roc, method = 'delong')

bmind_vs_enigma <- roc.test(bmind_roc, enigma_roc, method = 'delong')
bmind_vs_subtracted <- roc.test(bmind_roc, sub_roc, method = 'delong')

enigma_vs_subtracted <- roc.test(enigma_roc, sub_roc, method = 'delong')


