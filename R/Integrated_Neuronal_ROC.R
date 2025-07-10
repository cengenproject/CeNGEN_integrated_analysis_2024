### written 20240612 Alec Barrett
### updated 20240812 Alec Barrett
### updated again 20250706 Alec Barrett

### libraries

library(dplyr)
library(purrr)
library(ggplot2)
library(bayestestR)
library(pROC)
library(stringr)
library(patchwork)
library(pbapply)
library(pbmcapply)

### define functions ----
get_gt_metrics <- function(expression, ground_truth){
  
  
  expression <- expression |> as.matrix()
  expression_flat <- expression |> as.vector()
  ranker <- expression |> rank(ties.method = 'min') 
  orderer <- ranker |> order()
  helper <- ranker[orderer]
  thresholds <- unique(helper) |> as.numeric()
  
  
  ground_truth_ordered <- ground_truth |> as.matrix() |> as.vector()
  ground_truth_ordered <- ground_truth_ordered[orderer]
  
  t <- length(ground_truth_ordered)
  p <- sum(ground_truth_ordered)
  n <- t-p
  
  metrics <- sapply(thresholds, function(x){
    
    pos <- which(helper >= x)
    to <- length(pos)
    v <- ground_truth_ordered[pos]
    tp  <- sum(v)
    fp <- to - tp
    
    tpr <- tp/p
    fpr <- fp/n
    fdr <- fp/(tp+fp)
    
    return(c('TPR' = tpr,
             'FPR' = fpr,
             'FDR' = fdr))
  }) |> 
    t() |>
    as.data.frame()
  
  
  metrics$threshold <- expression_flat[orderer][thresholds]
  metrics
}

bootstrap_sensitivity <- function(n_boot, expr, gt, ncore) {
  
  pbmclapply(seq(1:n_boot), function(x){
    
    # Sample indices to maintain pairing between expr and gt
    set.seed(x)
    n_total <- length(gt |> as.matrix() |> as.vector())
    sample_indices <- sample(1:n_total, replace = T)
    
    gt_sample <- gt |> as.matrix() |> as.vector()
    gt_sample <- gt_sample[sample_indices]
    expr_sample <- expr |> as.matrix() |> as.vector() 
    expr_sample <- expr_sample[sample_indices]
    
    expression_flat <- expr_sample |> as.vector()
    ranker <- expression_flat |> rank(ties.method = 'min') 
    orderer <- ranker |> order()
    helper <- ranker[orderer]
    thresholds <- unique(helper) |> as.numeric()
    ground_truth_ordered <- gt_sample 
    ground_truth_ordered <- ground_truth_ordered[orderer]
    
    t <- length(ground_truth_ordered)
    p <- sum(ground_truth_ordered)
    n <- t-p
    
    diag <- sapply(thresholds, function(thresh){
      
      pos <- which(helper >= thresh)
      to <- length(pos)
      v <- ground_truth_ordered[pos]
      tp  <- sum(v)
      fp <- to - tp
      
      tpr <- tp/p
      fdr <- fp/(tp+fp)
      
      return(c('TPR' = tpr,
               'FDR' = fdr))
    }) |> 
      t() |>
      as.data.frame()
    
    diag <- rbind(diag, data.frame(TPR = 0, FDR = 0))
    
    if(0.05 %in% diag$FDR){
      returner <- diag[diag$FDR == 0.05,]
      returner <- returner[returner$TPR == max(returner$TPR),]
      if(sum(returner$TPR == max(returner$TPR)) > 1){
        returner <- returner[1,]
      }
      returner$predicted = 'no'
      
    }else{
      
      # Find closest points above and below 0.05
      above_05 <- diag[diag$FDR > 0.05,]
      below_05 <- diag[diag$FDR < 0.05,]
      
      if(nrow(above_05) == 0){
        # All FDR values are below 0.05, use the highest FDR
        returner <- diag[diag$FDR == max(diag$FDR),]
        returner$predicted = 'extrapolated'
      } else if(nrow(below_05) == 0){
        # All FDR values are above 0.05, use the lowest FDR  
        returner <- diag[diag$FDR == min(diag$FDR),]
        returner$predicted = 'extrapolated'
      } else {
        # Normal interpolation case
        point_above <- above_05[above_05$FDR == min(above_05$FDR),][1,]
        point_below <- below_05[below_05$FDR == max(below_05$FDR),][1,]
        
        interp_points <- rbind(point_above, point_below)
        TPR = approx(interp_points$FDR, interp_points$TPR, xout = 0.05)$y
        
        returner <- data.frame(TPR = TPR, FDR = 0.05, predicted = 'yes')
      }
    }
    return(returner)
    
  }, mc.cores = ncore)
}
#### load data


neuronal_gt <- read.table('references/bulk_all_ground_truth_121023.csv', sep = ',')
neuronal_gt$VD_DD <- neuronal_gt$VD

bulk_raw_TMM <- read.table('Data/bsn12_bulk_TMM_051624.tsv.gz', sep = '\t')
bulk_subtracted_TMM <- read.table('Data/bsn12_bulk_subtracted_TMM_070625.tsv.gz', sep = '\t')
bulk_integrated_aggregate <- read.table('Data/bsn12_subtracted_integrated_propadjust_070625.tsv')

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

## only consider common genes
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

## subset down for ground truth matrix
testing_gt <- neuronal_gt[neuronal_gt_genes, neurons]


## calculate Metrics across a wide range of thresholds

diags_aggr_raw_ave_plot <- tibble(get_gt_metrics(aggr_raw_TMM_plot, testing_gt),
                                  counts = "unaltered bulk")


diags_aggr_sub_TMM_plot <- tibble(get_gt_metrics(aggr_sub_TMM_plot, testing_gt),
                                  counts = "subtracted bulk")

diags_aggr_int_cpm_plot <- tibble(get_gt_metrics(bulk_integrated_aggregate_plot, testing_gt),
                                  counts = "integrated")

diags_proportions_plot <- tibble(get_gt_metrics(proportions_plot, testing_gt),
                                 counts = "sc proportions")

diags_sc_TPM_plot <- tibble(get_gt_metrics(sc_TPM_plot, testing_gt),
                            counts = "sc TPM")

diags_adjusted_proportions_plot <- tibble(get_gt_metrics(adjusted_proportions_plot, testing_gt),
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
ggsave('figures/Figure 4 Integrated analysis/B_Integrated_Neuronal_Testing_ROC_curves_070625.pdf', width = 9, height = 7)


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
ggsave('figures/Figure 5 Integrated analysis/B_Integrated_Neuronal_Testing_PR_curves_070625.pdf', width = 9, height = 7)



### threshold integrated data ----

##
integrated_spline <- splinefun(diags_aggr_int_cpm_plot$FDR,
                               diags_aggr_int_cpm_plot$threshold,
                               method = 'natural')

integrated_spline(x = 0.197)
integrated_spline(x = 0.14)
integrated_spline(x = 0.104)
integrated_spline(x = 0.084)


threshold_1_19.7p <- 0.1587915
threshold_2_14p <- 0.2824111
threshold_3_10.4p <- 0.4517556
threshold_4_8.4p <- 0.5805705


bulk_integrated_aggregate_threshold_1 <- bulk_integrated_aggregate
bulk_integrated_aggregate_threshold_2 <- bulk_integrated_aggregate
bulk_integrated_aggregate_threshold_3 <- bulk_integrated_aggregate
bulk_integrated_aggregate_threshold_4 <- bulk_integrated_aggregate

bulk_integrated_aggregate_threshold_1[bulk_integrated_aggregate_threshold_1 < threshold_1_19.7p] = 0
bulk_integrated_aggregate_threshold_2[bulk_integrated_aggregate_threshold_2 < threshold_2_14p] = 0
bulk_integrated_aggregate_threshold_3[bulk_integrated_aggregate_threshold_3 < threshold_3_10.4p] = 0
bulk_integrated_aggregate_threshold_4[bulk_integrated_aggregate_threshold_4 < threshold_4_8.4p] = 0


write.table(bulk_integrated_aggregate, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_unthresholded_070625.csv', 
            sep = ',',
            quote = F)
write.table(bulk_integrated_aggregate_threshold_1, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_1_070625.csv', 
            sep = ',',
            quote = F)
write.table(bulk_integrated_aggregate_threshold_2, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_2_070625.csv', 
            sep = ',',
            quote = F)
write.table(bulk_integrated_aggregate_threshold_3, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_3_070625.csv', 
            sep = ',',
            quote = F)
write.table(bulk_integrated_aggregate_threshold_4, 'Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_4_070625.csv', 
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
prop_tpm_delong$p.value * 10

roc_list_prop <- list('sc TPM' = sc_TPM_roc, 'sc proportions' = prop_roc)



single_cell_roc_df <- sapply(roc_list_prop, function(y){
  x <- ci.auc(y)
  return(c('lower_ci' = x[1], 'mean' = x[2], 'upper_ci' = x[3]))
}) |> t() |> data.frame() |> tibble::rownames_to_column('dataset') |>
  mutate(dataset = dataset |> factor(levels = c('sc TPM', 'sc proportions')))

ggplot(single_cell_roc_df) + 
  geom_col(aes(x = dataset, y = mean, fill = dataset), alpha = 0.8) + 
  geom_errorbar(aes(x = dataset, ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  #coord_cartesian(ylim = c(0.5,1)) +
  annotate("segment",
           x = 1,
           xend = 2,
           y = max(single_cell_roc_df$upper_ci) + 0.03,
           yend = max(single_cell_roc_df$upper_ci) + 0.03) +
  annotate("text", x = 1.5, y = max(single_cell_roc_df$upper_ci) + 0.05, 
           label = '*', 
           fontface = 'bold', size = 10) +
  geom_text(aes(x = dataset, y = mean+0.02, label = mean |> round(4)), fontface = 'bold') +
  ylab('AUROC') +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black', face = 'bold'), 
    axis.title = element_text(color = 'black', face = 'bold'),
    title = element_text(color = 'black', face = 'bold'))
ggsave('figures/Figure 4 Integrated analysis/C_Integrated_Neuronal_Testing_ROC_barchart_070625.pdf', width = 7, height = 7)



roc_list <- list('unaltered bulk' = raw_roc, 'subtracted bulk' = sub_roc, 'sc adjusted proportions' = adjusted_prop_roc,
                 'integrated' = int_roc)



ROC_df <- sapply(roc_list, function(y){
  x <- ci.auc(y)
  return(c('lower_ci' = x[1], 'mean' = x[2], 'upper_ci' = x[3]))
}) |> t() |> data.frame() |> tibble::rownames_to_column('dataset') |>
  mutate(dataset = dataset |> factor(levels = c('unaltered bulk', 'subtracted bulk', 'sc adjusted proportions', 'integrated')))


ggplot(ROC_df) + 
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
ggsave('figures/Figure 4 Integrated analysis/C_Integrated_Neuronal_Testing_ROC_barchart_070625.pdf', width = 7, height = 7)



## calculate the sensitivity of each dataset at 5% FDR, using bootstrap approach


system.time(unaltered_bulk_.05FDR_boot <- bootstrap_sensitivity(1000, aggr_raw_TMM_plot, testing_gt, 4))


subtracted_bulk_.05FDR_boot <- bootstrap_sensitivity(1000, aggr_sub_TMM_plot, testing_gt, 4)

integrated_.05FDR_boot <- bootstrap_sensitivity(1000, bulk_integrated_aggregate_plot, testing_gt, 4)

adjusted_proportions_.05FDR_boot <- bootstrap_sensitivity(1000, adjusted_proportions_plot, testing_gt, 4)



rbind(do.call(rbind, adjusted_proportions_.05FDR_boot) |> data.frame() |> mutate(dataset = 'sc adjusted proportions'),
      do.call(rbind, subtracted_bulk_.05FDR_boot) |> data.frame() |> mutate(dataset = 'subtracted bulk'),
      do.call(rbind, unaltered_bulk_.05FDR_boot) |> data.frame() |> mutate(dataset = 'unaltered bulk'),
      do.call(rbind, integrated_.05FDR_boot) |> data.frame() |> mutate(dataset = 'integrated')
) |> 
  data.frame() |> 
  mutate(dataset = factor(dataset, levels = c('unaltered bulk', 'subtracted bulk',
                                              'sc adjusted proportions', 'integrated'))) |>
  ggplot() +
  geom_boxplot(aes(x = dataset, y = TPR, fill = dataset), notch = T) +
  ggtitle('Sensitivity at 5% FDR') +
  xlab('') +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = 'black', face = 'bold'), 
        axis.title = element_text(color = 'black', face = 'bold'),
        title = element_text(color = 'black', face = 'bold'))
ggsave('figures/Figure 5 Integrated analysis/C_Sensitivity_at_5percent_FDR_070625.pdf', width = 7, height = 7)




