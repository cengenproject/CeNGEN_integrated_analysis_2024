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

bulk_raw_TMM <- read.table('Data/bsn12_bulk_TMM_051624.tsv.gz', sep = '\t')
bulk_subtracted_TMM <- read.table('Data/bsn12_bulk_subtracted_TMM_051624.tsv.gz', sep = '\t')
bulk_integrated_aggregate <- read.table('Data/bsn12_subtracted_integrated_propadjust_071724.tsv.gz')



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



## average bulk replicates per cell-type

aggr_raw_TMM <- bulk_raw_TMM
colnames(aggr_raw_TMM) <-str_split_fixed(colnames(aggr_raw_TMM),"r",2)[,1]
colnames(aggr_raw_TMM)[colnames(aggr_raw_TMM) %in% c('VD', 'DD')] <- 'VD_DD'
aggr_raw_TMM <- aggr_raw_TMM[,order(colnames(aggr_raw_TMM))]
aggr_raw_TMM <- data.frame(vapply(unique(colnames(aggr_raw_TMM)), function(x)
  rowMeans(aggr_raw_TMM[,colnames(aggr_raw_TMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_TMM)) ))
dim(aggr_raw_TMM)


aggr_subtracted_TMM <- bulk_subtracted_TMM
colnames(aggr_subtracted_TMM) <-str_split_fixed(colnames(aggr_subtracted_TMM),"r",2)[,1]
colnames(aggr_subtracted_TMM)[colnames(aggr_subtracted_TMM) %in% c('VD', 'DD')] <- 'VD_DD'
aggr_subtracted_TMM <- aggr_subtracted_TMM[,order(colnames(aggr_subtracted_TMM))]
aggr_subtracted_TMM <- data.frame(vapply(unique(colnames(aggr_subtracted_TMM)), function(x)
  rowMeans(aggr_subtracted_TMM[,colnames(aggr_subtracted_TMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_subtracted_TMM)) ))




## only consider common neurons

neurons <- intersect(colnames(aggr_subtracted_TMM), colnames(neuronal_gt))
neurons <- intersect(colnames(bulk_integrated_aggregate), neurons)


neuronal_gt_genes <- intersect(rownames(aggr_subtracted_TMM), rownames(neuronal_gt))
nrow(neuronal_gt)
sum(rownames(neuronal_gt) %in% rownames(aggr_raw_TMM))

### subset down to selected genes and neurons
aggr_raw_TMM_plot <- aggr_raw_TMM[neuronal_gt_genes, neurons]
aggr_sub_TMM_plot <- aggr_subtracted_TMM[neuronal_gt_genes, neurons]
sc_TPM_plot <- sc_TPM[neuronal_gt_genes,neurons]
proportions_plot <- prop_by_type[neuronal_gt_genes, neurons]
adjusted_proportions_plot <- prop_by_type_adjusted[neuronal_gt_genes, neurons]
bulk_integrated_aggregate_plot <- bulk_integrated_aggregate[neuronal_gt_genes, neurons]


testing_gt <- neuronal_gt[neuronal_gt_genes, neurons]


## calculate Metrics across thresholds

diags_aggr_raw_ave_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                                  TPR = map_dbl(threshold, ~get_tpr(aggr_raw_TMM_plot, testing_gt, .x)),
                                  FPR = map_dbl(threshold, ~get_fpr(aggr_raw_TMM_plot, testing_gt, .x)),
                                  FDR = map_dbl(threshold, ~get_fdr(aggr_raw_TMM_plot, testing_gt, .x)),
                                  counts = "unaltered bulk")


diags_aggr_sub_TMM_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                                  TPR = map_dbl(threshold, ~get_tpr(aggr_sub_TMM_plot, testing_gt, .x)),
                                  FPR = map_dbl(threshold, ~get_fpr(aggr_sub_TMM_plot, testing_gt, .x)),
                                  FDR = map_dbl(threshold, ~get_fdr(aggr_sub_TMM_plot, testing_gt, .x)),
                                  counts = "subtracted bulk")

diags_aggr_int_cpm_plot <- tibble(threshold = c(0,2**seq(-17,15,0.05)),
                                  TPR = map_dbl(threshold, ~get_tpr(bulk_integrated_aggregate_plot, testing_gt, .x)),
                                  FPR = map_dbl(threshold, ~get_fpr(bulk_integrated_aggregate_plot, testing_gt, .x)),
                                  FDR = map_dbl(threshold, ~get_fdr(bulk_integrated_aggregate_plot, testing_gt, .x)),
                                  counts = "integrated")

diags_proportions_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                                 TPR = map_dbl(threshold, ~get_tpr(proportions_plot, testing_gt, .x)),
                                 FPR = map_dbl(threshold, ~get_fpr(proportions_plot, testing_gt, .x)),
                                 FDR = map_dbl(threshold, ~get_fdr(proportions_plot, testing_gt, .x)),
                                 counts = "sc proportions")

diags_sc_TPM_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                            TPR = map_dbl(threshold, ~get_tpr(sc_TPM_plot, testing_gt, .x)),
                            FPR = map_dbl(threshold, ~get_fpr(sc_TPM_plot, testing_gt, .x)),
                            FDR = map_dbl(threshold, ~get_fdr(sc_TPM_plot, testing_gt, .x)),
                            counts = "sc TPM")

diags_adjusted_proportions_plot <- tibble(threshold = c(0,2**seq(-17,12,0.05)),
                                          TPR = map_dbl(threshold, ~get_tpr(adjusted_proportions_plot, testing_gt, .x)),
                                          FPR = map_dbl(threshold, ~get_fpr(adjusted_proportions_plot, testing_gt, .x)),
                                          FDR = map_dbl(threshold, ~get_fdr(adjusted_proportions_plot, testing_gt, .x)),
                                          counts = "sc adjusted proportions")


## plot ROC & PR curves


bind_rows(diags_aggr_raw_ave_plot,
          diags_aggr_sub_TMM_plot,
          diags_aggr_int_cpm_plot,
          diags_adjusted_proportions_plot
) |> 
  mutate(counts = factor(counts, levels = c("unaltered bulk", "subtracted bulk", "sc adjusted proportions", "integrated"))) |>
  ggplot(aes(x = FPR, y=TPR, color= counts)) +
  geom_abline(slope = 1, linetype='dashed') +
  geom_path(linewidth = 2, alpha = 0.7) +
  ggtitle('ROC for neuronal\ntesting genes') +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(color = 'black', face = 'bold'), 
        axis.title = element_text(color = 'black', face = 'bold'),
        title = element_text(color = 'black', face = 'bold'))
ggsave('figures/Figure 5 Integrated analysis/B_Integrated_Neuronal_Testing_ROC_curves_061224.pdf', width = 9, height = 7)


bind_rows(diags_aggr_raw_ave_plot,
          diags_aggr_sub_TMM_plot,
          diags_aggr_int_cpm_plot,
          diags_adjusted_proportions_plot
) |> 
  mutate(counts = factor(counts, levels = c("unaltered bulk", "subtracted bulk", "sc adjusted proportions", "integrated"))) |>
  ggplot(aes(x = 1-FDR, y=TPR, color= counts)) +
  #geom_vline(xintercept = min(1-diags_aggr_raw_ave_plot$FDR)) +
  geom_path(linewidth = 2, alpha = 0.7) +
  geom_path(data = data.frame(x = c(min(1-diags_aggr_raw_ave_plot$FDR), min(1-diags_aggr_raw_ave_plot$FDR)), y = c(0,1)), 
            aes(x=x,y=y),
            color = 'black',
            linetype = 'dashed', inherit.aes = F) +
  ggtitle('PR Curve for neuronal\ntesting genes') +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(color = 'black', face = 'bold'), 
        axis.title = element_text(color = 'black', face = 'bold'),
        title = element_text(color = 'black', face = 'bold'))
ggsave('figures/Figure 5 Integrated analysis/B_Integrated_Neuronal_Testing_PR_curves_061224.pdf', width = 9, height = 7)



### threshold integrated data ----

threshold_1_19.7p <- 0.12879
threshold_2_14p <- 0.23774
threshold_3_10.4p <- 0.38021
threshold_4_8.4p <- 0.53458


diags_aggr_int_cpm_plot


bulk_integrated_aggregate_threshold_1 <- bulk_integrated_aggregate
bulk_integrated_aggregate_threshold_2 <- bulk_integrated_aggregate
bulk_integrated_aggregate_threshold_3 <- bulk_integrated_aggregate
bulk_integrated_aggregate_threshold_4 <- bulk_integrated_aggregate

bulk_integrated_aggregate_threshold_1[bulk_integrated_aggregate_threshold_1 < threshold_1_19.7p] = 0
bulk_integrated_aggregate_threshold_2[bulk_integrated_aggregate_threshold_2 < threshold_2_14p] = 0
bulk_integrated_aggregate_threshold_3[bulk_integrated_aggregate_threshold_3 < threshold_3_10.4p] = 0
bulk_integrated_aggregate_threshold_4[bulk_integrated_aggregate_threshold_4 < threshold_4_8.4p] = 0


write.table(bulk_integrated_aggregate, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_unthresholded.csv', 
            sep = ',',
            quote = F)
write.table(bulk_integrated_aggregate_threshold_1, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_1.csv', 
            sep = ',',
            quote = F)
write.table(bulk_integrated_aggregate_threshold_2, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_2.csv', 
            sep = ',',
            quote = F)
write.table(bulk_integrated_aggregate_threshold_3, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_3.csv', 
            sep = ',',
            quote = F)
write.table(bulk_integrated_aggregate_threshold_4, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_4.csv', 
            sep = ',',
            quote = F)




## perform statistical tests with DeLong test ----

aggr_raw_TMM_plot |> unlist() |> length()

raw_roc <- roc(testing_gt |> unlist(),
               aggr_raw_TMM_plot|> log1p() |> unlist())

sub_roc <- roc(testing_gt |> unlist(),
               aggr_sub_TMM_plot |> log1p() |> unlist())

int_roc <- roc(testing_gt |> unlist(),
               bulk_integrated_aggregate_plot |> log1p() |> unlist())

prop_roc <- roc(testing_gt |> unlist(),
                proportions_plot |> unlist())

sc_TPM_roc <- roc(testing_gt |> unlist(),
                  sc_TPM_plot |> unlist())

adjusted_prop_roc <- roc(testing_gt |> unlist(),
                         adjusted_proportions_plot |> unlist())

roc.test(raw_roc, sub_roc)
roc.test(raw_roc, prop_roc)
roc.test(raw_roc, adjusted_prop_roc)
roc.test(raw_roc, int_roc)

roc.test(sub_roc, prop_roc)
roc.test(sub_roc, adjusted_prop_roc)
roc.test(sub_roc, int_roc)

roc.test(prop_roc, adjusted_prop_roc)
roc.test(prop_roc, int_roc)

roc.test(adjusted_prop_roc, int_roc)

## plot barcharts of AUROC with CI


prop_tpm_delong <- roc.test(prop_roc, sc_TPM_roc)
prop_tpm_delong$p.value * 12

roc_list_prop <- list('sc TPM' = sc_TPM_roc, 'sc proportions' = prop_roc)



sapply(roc_list_prop, function(y){
  x <- ci.auc(y)
  return(c('lower_ci' = x[1], 'mean' = x[2], 'upper_ci' = x[3]))
}) |> t() |> data.frame() |> tibble::rownames_to_column('dataset') |>
  mutate(dataset = dataset |> factor(levels = c('sc TPM', 'sc proportions'))) |>
  ggplot() + 
  geom_col(aes(x = dataset, y = mean, fill = dataset), alpha = 0.8) + 
  geom_errorbar(aes(x = dataset, ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  #coord_cartesian(ylim = c(0.5,1)) +
  geom_text(aes(x = dataset, y = mean+0.02, label = mean |> round(4)), fontface = 'bold') +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black', face = 'bold'), 
    axis.title = element_text(color = 'black', face = 'bold'),
    title = element_text(color = 'black', face = 'bold'))
ggsave('figures/Figure 4 Integrated analysis/C_Integrated_Neuronal_Testing_ROC_barchart_061224.pdf', width = 7, height = 7)



roc_list <- list('unaltered bulk' = raw_roc, 'subtracted bulk' = sub_roc, 'sc adjusted proportions' = adjusted_prop_roc,
                 'integrated' = int_roc)



sapply(roc_list, function(y){
  x <- ci.auc(y)
  return(c('lower_ci' = x[1], 'mean' = x[2], 'upper_ci' = x[3]))
}) |> t() |> data.frame() |> tibble::rownames_to_column('dataset') |>
  mutate(dataset = dataset |> factor(levels = c('unaltered bulk', 'subtracted bulk', 'sc adjusted proportions', 'integrated'))) |>
  ggplot() + 
  geom_col(aes(x = dataset, y = mean, fill = dataset), alpha = 0.8) + 
  geom_errorbar(aes(x = dataset, ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  #coord_cartesian(ylim = c(0.5,1)) +
  geom_text(aes(x = dataset, y = mean+0.02, label = mean |> round(4)), fontface = 'bold') +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black', face = 'bold'), 
    axis.title = element_text(color = 'black', face = 'bold'),
    title = element_text(color = 'black', face = 'bold'))
ggsave('figures/Figure 4 Integrated analysis/C_Integrated_Neuronal_Testing_ROC_barchart_061224.pdf', width = 7, height = 7)



## calculate the sensitivity of each dataset at 5% FDR, using bootstrap approach

unaltered_bulk_.05FDR_boot <- pblapply(seq(1:1000), function(x){
  
  
  ## sample with replacement
  set.seed(x)
  gt <- testing_gt |> unlist() |> sample(replace = T)
  
  set.seed(x)
  bin <- aggr_raw_TMM_plot |> unlist() |> sample(replace = T) |> log1p()
  
  
  ## caculate TPR/FDR across thresholds
  diag <- tibble(threshold = seq(5,10 ,0.001),
                 TPR = map_dbl(threshold, ~get_tpr(bin, gt, .x)),
                 FDR = map_dbl(threshold, ~get_fdr(bin, gt, .x))) |> data.frame()
  #diag <- unique(diag[,c('TPR', 'FDR')])
  
  
  ## if the dataset can hit exactly 5% FDR, then use that value, but if not, make a linear model between the two surrounding points only, and interpolate the TPR at 5% FDR
  if(0.05 %in% diag$FDR){
    returner <- diag[diag$FDR == 0.05,]
    returner <- returner[returner$TPR == max(returner$TPR),]
    if(sum(returner$TPR == max(returner$TPR)) > 1){
      returner <- returner[1,]
    }
    returner$predicted = 'no'
    
  }else{
    
    diag$distance <- 0.05-diag$FDR
    
    diag_above <- diag[diag$distance < 0,]
    diag_above <- diag_above[diag_above$distance == max(diag_above$distance),]
    diag_below <- diag[diag$distance > 0,]
    if(sum(diag$distance > 0)==0){
      diag_below <- data.frame(threshold = Inf, TPR = 0, FDR = 0, distance = 0.05)
    }
    else{diag_below <- diag_below[diag_below$distance == min(diag_below$distance),]}
    diag <- rbind(diag_above, diag_below)
    
    FDR = 0.05
    TPR = predict(lm(TPR ~ FDR, diag), list('FDR' = 0.05)) [[1]]
    
    returner <- data.frame(threshold = mean(diag$threshold), TPR = TPR, FDR = FDR, predicted = 'yes')
    
  }
  return(returner)
  
})

integrated_.05FDR_boot <- pblapply(seq(1:1000), function(x){
  
  set.seed(x)
  gt <- testing_gt |> unlist() |> sample(replace = T)
  
  set.seed(x)
  bin <- bulk_integrated_aggregate_plot |> unlist() |> sample(replace = T) |> log1p()
  
  diag <- tibble(threshold = seq(0,2 ,0.001),
                 TPR = map_dbl(threshold, ~get_tpr(bin, gt, .x)),
                 FDR = map_dbl(threshold, ~get_fdr(bin, gt, .x))) |> data.frame()
  #diag <- unique(diag[,c('TPR', 'FDR')])
  
  
  ## if the dataset can hit exactly 5% FDR, then use that value, but if not, make a linear model between the two surrounding points only, and interpolate the TPR at 5% FDR
  if(0.05 %in% diag$FDR){
    returner <- diag[diag$FDR == 0.05,]
    returner <- returner[returner$TPR == max(returner$TPR),]
    if(sum(returner$TPR == max(returner$TPR)) > 1){
      returner <- returner[1,]
    }
    returner$predicted = 'no'
    
  }else{
    
    diag$distance <- 0.05-diag$FDR
    
    diag_above <- diag[diag$distance < 0,]
    diag_above <- diag_above[diag_above$distance == max(diag_above$distance),]
    diag_below <- diag[diag$distance > 0,]
    if(sum(diag$distance > 0)==0){
      diag_below <- data.frame(threshold = Inf, TPR = 0, FDR = 0, distance = 0.05)
    }
    else{diag_below <- diag_below[diag_below$distance == min(diag_below$distance),]}
    diag <- rbind(diag_above, diag_below)
    
    FDR = 0.05
    TPR = predict(lm(TPR ~ FDR, diag), list('FDR' = 0.05)) [[1]]
    
    returner <- data.frame(threshold = mean(diag$threshold), TPR = TPR, FDR = FDR, predicted = 'yes')
    
  }
  return(returner)
  
})


subtracted_bulk_.05FDR_boot <- pblapply(seq(1:1000), function(x){
  
  set.seed(x)
  gt <- testing_gt |> unlist() |> sample(replace = T)
  
  set.seed(x)
  bin <- aggr_sub_TMM_plot |> unlist() |> sample(replace = T) |> log1p()
  
  diag <- tibble(threshold = seq(4,6 ,0.001),
                 TPR = map_dbl(threshold, ~get_tpr(bin, gt, .x)),
                 FDR = map_dbl(threshold, ~get_fdr(bin, gt, .x))) |> data.frame()
  #diag <- unique(diag[,c('TPR', 'FDR')])
  
  
  ## if the dataset can hit exactly 5% FDR, then use that value, but if not, make a linear model between the two surrounding points only, and interpolate the TPR at 5% FDR
  if(0.05 %in% diag$FDR){
    returner <- diag[diag$FDR == 0.05,]
    returner <- returner[returner$TPR == max(returner$TPR),]
    if(sum(returner$TPR == max(returner$TPR)) > 1){
      returner <- returner[1,]
    }
    returner$predicted = 'no'
    
  }else{
    
    diag$distance <- 0.05-diag$FDR
    
    diag_above <- diag[diag$distance < 0,]
    diag_above <- diag_above[diag_above$distance == max(diag_above$distance),]
    diag_below <- diag[diag$distance > 0,]
    if(sum(diag$distance > 0)==0){
      diag_below <- data.frame(threshold = Inf, TPR = 0, FDR = 0, distance = 0.05)
    }
    else{diag_below <- diag_below[diag_below$distance == min(diag_below$distance),]}
    diag <- rbind(diag_above, diag_below)
    
    FDR = 0.05
    TPR = predict(lm(TPR ~ FDR, diag), list('FDR' = 0.05)) [[1]]
    
    returner <- data.frame(threshold = mean(diag$threshold), TPR = TPR, FDR = FDR, predicted = 'yes')
    
  }
  return(returner)
  
})
subtracted_bulk_.05FDR_boot




adjusted_proportions_.05FDR_boot <- pblapply(seq(1:1000), function(x){
  
  set.seed(x)
  gt <- testing_gt |> unlist() |> sample(replace = T)
  
  set.seed(x)
  bin <- adjusted_proportions_plot |> unlist() |> sample(replace = T) 
  
  diag <- tibble(threshold = seq(0.1,0.5 ,0.0005),
                 TPR = map_dbl(threshold, ~get_tpr(bin, gt, .x)),
                 FDR = map_dbl(threshold, ~get_fdr(bin, gt, .x))) |> data.frame()
  #diag <- unique(diag[,c('TPR', 'FDR')])
  
  
  
  ## if the dataset can hit exactly 5% FDR, then use that value, but if not, make a linear model between the two surrounding points only, and interpolate the TPR at 5% FDR
  if(0.05 %in% diag$FDR){
    returner <- diag[diag$FDR == 0.05,]
    returner <- returner[returner$TPR == max(returner$TPR),]
    if(sum(returner$TPR == max(returner$TPR)) > 1){
      returner <- returner[1,]
    }
    returner$predicted = 'no'
    
  }else{
    
    diag$distance <- 0.05-diag$FDR
    
    diag_above <- diag[diag$distance < 0,]
    diag_above <- diag_above[diag_above$distance == max(diag_above$distance),]
    diag_below <- diag[diag$distance > 0,]
    if(sum(diag$distance > 0)==0){
      diag_below <- data.frame(threshold = Inf, TPR = 0, FDR = 0, distance = 0.05)
    }
    else{diag_below <- diag_below[diag_below$distance == min(diag_below$distance),]}
    diag <- rbind(diag_above, diag_below)
    
    FDR = 0.05
    TPR = predict(lm(TPR ~ FDR, diag), list('FDR' = 0.05)) [[1]]
    
    returner <- data.frame(threshold = mean(diag$threshold), TPR = TPR, FDR = FDR, predicted = 'yes')
    
  }
  return(returner)
  
})

rbind(do.call(rbind, adjusted_proportions_.05FDR_boot) |> data.frame() |> mutate(dataset = 'sc adjusted proportions'),
      do.call(rbind, subtracted_bulk_.05FDR_boot) |> data.frame() |> mutate(dataset = 'subtracted bulk'),
      do.call(rbind, unaltered_bulk_.05FDR_boot) |> data.frame() |> mutate(dataset = 'unaltered bulk'),
      do.call(rbind, integrated_.05FDR_boot) |> data.frame() |> mutate(dataset = 'integrated')
) |> data.frame() |> mutate(dataset = factor(dataset, levels = c('unaltered bulk', 'subtracted bulk', 'sc adjusted proportions', 'integrated'))) |>
  ggplot() +
  geom_boxplot(aes(x = dataset, y = TPR, fill = dataset), notch = T) +
  ggtitle('Sensitivity at 5% FDR') +
  xlab('') +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = 'black', face = 'bold'), 
        axis.title = element_text(color = 'black', face = 'bold'),
        title = element_text(color = 'black', face = 'bold'))
ggsave('figures/Figure 5 Integrated analysis/C_Sensitivity_at_5percent_FDR_240614.pdf', width = 7, height = 7)




