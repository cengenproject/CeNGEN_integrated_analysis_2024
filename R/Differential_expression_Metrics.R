### check qlf output -- expected to run after "single_cell_differential_expression.R"

library(dplyr)
library(wbData)
library(edgeR)
library(pbapply)
library(reshape)
library(stringr)
library(ggplot2)
library(reshape)


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

get_gt_statistics <- function(qlf_table, gt_df,
                              cell1, cell2,
                              logFC_threshold = 1,
                              pvalue_name = 'p.adj', pvalue_threshold = 0.05){
  
  gt_df <- gt[,c(cell1, cell2)]
  gt_df <- gt_df[rowSums(gt_df == 1) < 2,]
  
  truth_sign <- sign(gt_df[,cell1] - gt_df[,cell2])
  names(truth_sign) <- rownames(gt_df)
  
  rn <- sum(truth_sign == 0)
  rp <- sum(truth_sign != 0)
  
  qlf_table <- qlf_table[rownames(gt_df),]
  rownames(qlf_table) <- rownames(gt_df)
  qlf_table[is.na(qlf_table$logFC),] <- 0.1
  
  qlf_table$logFC_thresholded <- qlf_table$logFC
  qlf_table$logFC_thresholded[abs(qlf_table$logFC_thresholded) < logFC_threshold] <- 0
  
  
  
  tp <- sum((truth_sign == sign(qlf_table$logFC_thresholded) & 
               qlf_table[[pvalue_name]] < pvalue_threshold) &
              truth_sign != 0 )
  tn <- sum((qlf_table$logFC_thresholded == 0 |
               qlf_table[[pvalue_name]] > pvalue_threshold) &
              truth_sign == 0)
  
  fp <- sum((truth_sign != sign(qlf_table$logFC_thresholded)) & qlf_table$logFC_thresholded != 0 &
              qlf_table[[pvalue_name]] < pvalue_threshold)
  
  fn <- sum((qlf_table$logFC_thresholded == 0 |
               qlf_table[[pvalue_name]] > pvalue_threshold) & truth_sign != 0)
  
  return(c('TPR' = tp/rp,
           'TNR' = tn/rn,
           'FPR' = fp/rn,
           'FNR' = fn/rp,
           'FDR' = fp/(tp+fp),
           'TDR' = tp/(tp+fp),
           'MCC' = (((tp*tn) - (fp*fn))) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))))
  
  
}





## load data

gt <- read.table('references/bulk_all_ground_truth_121023.csv', sep = ',')
gt$VD_DD <- gt$VD

sc_counts_qlf_list <- readRDS('Data_out/sc_counts_qlf_tables.rds')
sc_p2c_qlf_list <- readRDS('Data_out/sc_prop_qlf_tables.rds')

### compare one 



sc_counts_qlf_list <- pblapply(sc_counts_qlf_list, function(x){
  
  x$FDR <- p.adjust(x$PValue, 'BH')
  x$FWER <- p.adjust(x$PValue, 'bonferroni')
  
  return(x)
})

sc_p2c_qlf_list <- pblapply(sc_p2c_qlf_list, function(x){
  
  x$FDR <- p.adjust(x$PValue, 'BH')
  x$FWER <- p.adjust(x$PValue, 'bonferroni')
  
  return(x)
})



sc_counts_stats <- pbsapply(names(sc_counts_qlf_list), function(x){
  
  cell1 <- str_split_fixed(gsub('\\(|\\)', '', x), '-', 2)[,1]
  cell2 <- str_split_fixed(gsub('\\(|\\)', '', x), '-', 2)[,2]
  
  return(get_gt_statistics(sc_counts_qlf_list[[x]], gt, cell1, cell2, 
                           logFC_threshold = 0.5, pvalue_name = 'FWER'))
  
  
}) |> t() |> data.frame(dataset = 'sc counts')

sc_p2c_stats <- pbsapply(names(sc_p2c_qlf_list), function(x){
  
  cell1 <- str_split_fixed(gsub('\\(|\\)', '', x), '-', 2)[,1]
  cell2 <- str_split_fixed(gsub('\\(|\\)', '', x), '-', 2)[,2]
  
  return(get_gt_statistics(sc_p2c_qlf_list[[x]], gt, cell1, cell2, 
                           logFC_threshold = 0.5, pvalue_name = 'FWER'))
  
  
}) |> t() |> data.frame(dataset = 'p2c counts')


rbind(sc_counts_stats |> melt(), sc_p2c_stats |> melt()) |> 
  filter(variable %in% c('FDR', 'TPR', 'FPR', 'MCC')) |>
  mutate(dataset = factor(dataset, levels = c('sc ranksum', 'sc counts', 'p2c counts'))) |>
  ggplot() +
  geom_boxplot(aes(x = variable, y = value, fill = dataset), alpha = 0.7, notch = T)
ggsave('figures/Differential_Expression_metrics.pdf', width = 7, heigh = 5)


wilcox.test(sc_counts_stats$TPR, sc_p2c_stats$TPR, paired = F)$p.value * 4
wilcox.test(sc_counts_stats$FPR, sc_p2c_stats$FPR, paired = F)$p.value * 4
wilcox.test(sc_counts_stats$FDR, sc_p2c_stats$FDR, paired = F)$p.value * 4
wilcox.test(sc_counts_stats$MCC, sc_p2c_stats$MCC, paired = F)$p.value * 4


