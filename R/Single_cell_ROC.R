### libraries

library(dplyr)
library(purrr)
library(ggplot2)
library(bayestestR)
library(pROC)
library(stringr)
library(patchwork)



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


neuronal_gt <- read.table('references/bulk_all_ground_truth_121023.csv', sep = ',')
neuronal_gt$VD_DD <- neuronal_gt$VD


sc_TPM <- read.table('Data/CeNGEN_TPM_080421.tsv.gz')
sc_TPM$VD_DD <- sc_TPM$VD


prop_by_type <- read.table('Data/SingleCell_proportions_Bulk_annotations.tsv.gz')
colnames(prop_by_type)
prop_by_type <- prop_by_type[,setdiff(colnames(prop_by_type), c('Intestine',
                                                                'Glia',
                                                                'Pharynx',
                                                                'Excretory',
                                                                'Hypodermis',
                                                                'Rectal_cells',
                                                                'Reproductive',
                                                                'Muscle_mesoderm'))]



# hard thresholds fixed
low_hard <- 0.02
high_hard <- 0.01


prop_by_type_adjusted <- prop_by_type

prop_by_type_adjusted[apply(prop_by_type_adjusted, 1, min) > high_hard, ] <- 1
prop_by_type_adjusted[apply(prop_by_type_adjusted, 1, max) < low_hard, ] <- 0
mid_range <- (apply(prop_by_type_adjusted, 1, max) > low_hard) & (apply(prop_by_type_adjusted, 1, min) < high_hard)
prop_by_type_adjusted[mid_range, ] <-
  prop_by_type_adjusted[mid_range, ]/apply(prop_by_type_adjusted[mid_range,], 1, max)



neurons <- c("ADL", "AFD", "AIM", "AIN", "AIY", "ASEL", "ASER", 
             "ASG", "ASI", "ASK", "AVA", "AVE", "AVG", "AVH", "AVK",
             "AVL", "AVM", "AWA", "AWB", "AWC", "BAG", "CAN", "CEP", 
             "DA", "DVB", "DVC", "HSN", "I5", "IL1", "IL2", "LUA", "NSM",
             "OLL", "OLQ", "PHA", "PVC", "PVD", "PVM", "PVP", "PVQ",
             "RIA", "RIC", "RIM", "RIS", "RMD", "RME", "SIA", "SMB",
             "SMD", "VB", "VC", "VD_DD")

neuronal_gt_genes <- rownames(neuronal_gt)


sc_TPM_plot <- sc_TPM[neuronal_gt_genes,neurons]
proportions_plot <- prop_by_type[neuronal_gt_genes, neurons]

testing_gt <- neuronal_gt[neuronal_gt_genes, neurons]


diags_sc_TPM_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                            TPR = map_dbl(threshold, ~get_tpr(sc_TPM_plot, testing_gt, .x)),
                            FPR = map_dbl(threshold, ~get_fpr(sc_TPM_plot, testing_gt, .x)),
                            FDR = map_dbl(threshold, ~get_fdr(sc_TPM_plot, testing_gt, .x)),
                            counts = "sc TPM")

diags_proportions_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                                 TPR = map_dbl(threshold, ~get_tpr(proportions_plot, testing_gt, .x)),
                                 FPR = map_dbl(threshold, ~get_fpr(proportions_plot, testing_gt, .x)),
                                 FDR = map_dbl(threshold, ~get_fdr(proportions_plot, testing_gt, .x)),
                                 counts = "sc proportions")



## plot ROC & PR curves

bind_rows(diags_sc_TPM_plot,
          diags_proportions_plot) |> 
  mutate(counts = factor(counts, levels = c("sc TPM", 'sc proportions', "sc adjusted proportions"))) |>
  ggplot(aes(x = FPR, y=TPR, color= counts)) +
  geom_abline(slope = 1, linetype='dashed') +
  geom_path(linewidth = 2, alpha = 0.7) +
  ggtitle('ROC for neuronal\ntesting genes') +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(color = 'black', face = 'bold'), 
        axis.title = element_text(color = 'black', face = 'bold'),
        title = element_text(color = 'black', face = 'bold'))
ggsave('figures/Prop2Count/Ai_proportions_vs_TPM_ROC_240616.pdf', width = 7, height = 5)
