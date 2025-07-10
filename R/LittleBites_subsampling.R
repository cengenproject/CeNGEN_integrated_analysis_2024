# Load required libraries
library(Seurat)
library(dplyr)
library(caret)
library(LittleBites)
library(nnls)
library(parallel)
library(edgeR)
library(pbmcapply)
library(broom)



###functions 

calculate_cv_results <- function(deconv_results, corr_results) {
  cv_results <- list()
  
  for(sample_name in names(deconv_results)) {
    cat(paste("Calculating CV for sample:", sample_name, "\n"))
    
    cv_results[[sample_name]] <- list()
    
    for(subsample_size in names(deconv_results[[sample_name]])) {
      
      # Extract metrics across all iterations
      initial_target_purity <- sapply(deconv_results[[sample_name]][[subsample_size]],
                                      function(x) x$initial_estimates[1])  # First element is target
      final_target_purity <- sapply(deconv_results[[sample_name]][[subsample_size]],
                                    function(x) x$final_estimates[1])
      train_aurocs <- sapply(deconv_results[[sample_name]][[subsample_size]],
                             function(x) x$train_AUROC)
      test_aurocs <- sapply(deconv_results[[sample_name]][[subsample_size]],
                            function(x) x$test_AUROC)
      
      # Extract subtracted value statistics
      subtracted_means <- sapply(deconv_results[[sample_name]][[subsample_size]],
                                 function(x) x$subtracted_mean)
      subtracted_medians <- sapply(deconv_results[[sample_name]][[subsample_size]],
                                   function(x) x$subtracted_median)
      
      # NEW: Extract correlation and gene stability metrics
      mean_correlation <- corr_results[[sample_name]][[subsample_size]]$correlation_metrics$mean_pairwise_correlation
      n_correlation_pairs <- corr_results[[sample_name]][[subsample_size]]$correlation_metrics$n_pairs
      
      # Gene stability metrics
      mean_gene_cv <- corr_results[[sample_name]][[subsample_size]]$gene_stability_metrics$mean_gene_cv
      median_gene_cv <- corr_results[[sample_name]][[subsample_size]]$gene_stability_metrics$median_gene_cv
      stable_gene_proportion <- corr_results[[sample_name]][[subsample_size]]$gene_stability_metrics$stable_gene_proportion
      cv_90th_percentile <- corr_results[[sample_name]][[subsample_size]]$gene_stability_metrics$cv_90th_percentile
      
      # Calculate coefficient of variation (CV = sd/mean)
      cv_results[[sample_name]][[subsample_size]] <- list(
        initial_purity_mean = mean(initial_target_purity, na.rm = TRUE),
        initial_purity_cv = sd(initial_target_purity, na.rm = TRUE) / mean(initial_target_purity, na.rm = TRUE),
        final_purity_mean = mean(final_target_purity, na.rm = TRUE),
        final_purity_cv = sd(final_target_purity, na.rm = TRUE) / mean(final_target_purity, na.rm = TRUE),
        train_auroc_mean = mean(train_aurocs, na.rm = TRUE),
        train_auroc_cv = sd(train_aurocs, na.rm = TRUE) / mean(train_aurocs, na.rm = TRUE),
        test_auroc_mean = mean(test_aurocs, na.rm = TRUE),
        test_auroc_cv = sd(test_aurocs, na.rm = TRUE) / mean(test_aurocs, na.rm = TRUE),
        # Subtracted value statistics
        subtracted_mean_mean = mean(subtracted_means, na.rm = TRUE),
        subtracted_mean_cv = sd(subtracted_means, na.rm = TRUE) / mean(subtracted_means, na.rm = TRUE),
        subtracted_median_mean = mean(subtracted_medians, na.rm = TRUE),
        subtracted_median_cv = sd(subtracted_medians, na.rm = TRUE) / mean(subtracted_medians, na.rm = TRUE),
        # NEW: Correlation metrics
        mean_pairwise_correlation = mean_correlation,
        n_correlation_pairs = n_correlation_pairs,
        # NEW: Gene stability metrics
        mean_gene_cv = mean_gene_cv,
        median_gene_cv = median_gene_cv,
        stable_gene_proportion = stable_gene_proportion,
        cv_90th_percentile = cv_90th_percentile,
        n_iterations = length(initial_target_purity)
      )
    }
  }
  
  return(cv_results)
}

calculate_pairwise_correlations <- function(subtracted_matrix) {
  
  log_transformed <- log1p(subtracted_matrix)
  correlation_matrix <- cor(log_transformed, use = "complete.obs")
  
  # Extract upper triangle (excluding diagonal) to get unique pairwise correlations
  upper_triangle <- upper.tri(correlation_matrix, diag = FALSE)
  pairwise_correlations <- correlation_matrix[upper_triangle]
  
  # Calculate mean pairwise correlation
  mean_pairwise_correlation <- mean(pairwise_correlations, na.rm = TRUE)
  
  return(list(
    pairwise_correlations = pairwise_correlations,
    mean_pairwise_correlation = mean_pairwise_correlation,
    n_pairs = length(pairwise_correlations)
  ))
}

calculate_gene_stability <- function(subtracted_matrix) {
  gene_cv <- apply(subtracted_matrix, 1, function(x) {
    gene_mean <- mean(x, na.rm=TRUE)
    if(gene_mean == 0) return(NA)  # Handle zero-mean genes
    sd(x, na.rm=TRUE) / gene_mean
  })
  
  return(list(
    mean_gene_cv = mean(gene_cv, na.rm=TRUE),
    median_gene_cv = median(gene_cv, na.rm=TRUE),
    stable_gene_proportion = sum(gene_cv < 0.2, na.rm=TRUE) / sum(!is.na(gene_cv)),
    cv_90th_percentile = quantile(gene_cv, 0.9, na.rm=TRUE)
  ))
}

create_cv_summary_table <- function(cv_results) {
  summary_df <- data.frame()
  
  for(sample_name in names(cv_results)) {
    for(subsample_size in names(cv_results[[sample_name]])) {
      
      row_data <- data.frame(
        sample = sample_name,
        neuron_type = str_split_fixed(sample_name, 'r', 2)[,1],
        subsample_size = as.numeric(subsample_size),
        initial_purity_mean = cv_results[[sample_name]][[subsample_size]]$initial_purity_mean,
        initial_purity_cv = cv_results[[sample_name]][[subsample_size]]$initial_purity_cv,
        final_purity_mean = cv_results[[sample_name]][[subsample_size]]$final_purity_mean,
        final_purity_cv = cv_results[[sample_name]][[subsample_size]]$final_purity_cv,
        train_auroc_mean = cv_results[[sample_name]][[subsample_size]]$train_auroc_mean,
        train_auroc_cv = cv_results[[sample_name]][[subsample_size]]$train_auroc_cv,
        test_auroc_mean = cv_results[[sample_name]][[subsample_size]]$test_auroc_mean,
        test_auroc_cv = cv_results[[sample_name]][[subsample_size]]$test_auroc_cv,
        # Subtracted value statistics
        subtracted_mean_mean = cv_results[[sample_name]][[subsample_size]]$subtracted_mean_mean,
        subtracted_mean_cv = cv_results[[sample_name]][[subsample_size]]$subtracted_mean_cv,
        subtracted_median_mean = cv_results[[sample_name]][[subsample_size]]$subtracted_median_mean,
        subtracted_median_cv = cv_results[[sample_name]][[subsample_size]]$subtracted_median_cv,
        # NEW: Correlation metrics
        mean_pairwise_correlation = cv_results[[sample_name]][[subsample_size]]$mean_pairwise_correlation,
        n_correlation_pairs = cv_results[[sample_name]][[subsample_size]]$n_correlation_pairs,
        # NEW: Gene stability metrics
        mean_gene_cv = cv_results[[sample_name]][[subsample_size]]$mean_gene_cv,
        median_gene_cv = cv_results[[sample_name]][[subsample_size]]$median_gene_cv,
        stable_gene_proportion = cv_results[[sample_name]][[subsample_size]]$stable_gene_proportion,
        cv_90th_percentile = cv_results[[sample_name]][[subsample_size]]$cv_90th_percentile,
        stringsAsFactors = FALSE
      )
      
      summary_df <- rbind(summary_df, row_data)
    }
  }
  
  return(summary_df)
}


#### littlebites reference degradation testing

sc_ob <- readRDS('../single_cell_data/032224_L4_all_cells_Seurat5.rds')
### seurat object downloaded from CeNGEN website downloads page, ~0.5 gb disk space, ~ 2gb RAM to open

sc_ob_counts <- sc_ob@assays$RNA@counts
sc_ob_meta <- sc_ob@meta.data

sc_ob <- CreateSeuratObject(counts = sc_ob_counts, meta.data = sc_ob_meta)

sc_ob <- sc_ob[,sc_ob$Tissue != 'Unknown' & sc_ob$Tissue != 'Unannotated']


sc_ob@meta.data$Cell.type[sc_ob@meta.data$Cell.type %in% c('DB01')] <- 'DB'
sc_ob@meta.data$Cell.type[sc_ob@meta.data$Cell.type %in% c('DA9')] <- 'DA'
sc_ob@meta.data$Cell.type[sc_ob@meta.data$Cell.type %in% c('VC_4_5')] <- 'VC'
sc_ob@meta.data$Cell.type[sc_ob@meta.data$Cell.type %in% c('VB01', 'VB02')] <- 'VB'
sc_ob@meta.data$Cell.type[sc_ob@meta.data$Cell.type %in% c('RMD_DV', 'RMD_LR')] <- 'RMD'
sc_ob@meta.data$Cell.type[sc_ob@meta.data$Cell.type %in% c('RME_DV', 'RME_LR')] <- 'RME'
sc_ob@meta.data$Cell.type[sc_ob@meta.data$Cell.type %in% c('VA12')] <- 'VA'
sc_ob@meta.data$Cell.type[sc_ob@meta.data$Cell.type %in% c('IL2_DV', 'IL2_LR')] <- 'IL2'
sc_ob@meta.data$Cell.type[sc_ob@meta.data$Cell.type %in% c('AWC_ON', 'AWC_OFF')] <- 'AWC'
sc_ob <- sc_ob[,!(sc_ob$Cell.type %in% c('RIV_stressed', 'SMD_stressed'))]
sc_ob$Tissue |> table()

non_neuronal_list <- c('Excretory', 'Glia', 'Hypodermis', 'Intestine', 'Muscle_mesoderm', 'Pharynx', 'Reproductive')


sc_ob$neuron_level <- sc_ob$Tissue

sc_ob$neuron_level[sc_ob$neuron_level=='Neuron'] <- sc_ob$Cell.type[sc_ob$neuron_level=='Neuron']


sc_ob$neuron_level |> table()


cells_to_keep <- c(non_neuronal_list,
                   'CAN', 'DB', 'VB', 'AVJ', 'RIG',
                   'RMD', 'SIA', 'SMB', 'SMD', 'VD_DD')

keep_mask <- sc_ob$neuron_level %in% cells_to_keep

sc_ob <- sc_ob[,keep_mask]


sc_ob


### ground truth ----
UNN_ground_truth_mtx <- read.table('references/ubituiqtous_and_nonNeuronal_gt_genes_matrix_042222.tsv', sep = '\t')
train_test_split <- readRDS('references/ubituiqtous_and_nonNeuronal_gt_genes_split_042222.rds')
UNN_train <- UNN_ground_truth_mtx[train_test_split$training_genes,]
UNN_train$DD <- UNN_train$VD_DD
UNN_train$VD <- UNN_train$VD_DD
UNN_test <- UNN_ground_truth_mtx[train_test_split$testing_genes,]
UNN_test$DD <- UNN_test$VD_DD
UNN_test$VD <- UNN_test$VD_DD

non_neuronal_types <- c("Glia", "Hypodermis", "Intestine",
                        "Muscle_mesoderm", "Pharynx", "Reproductive")

## get indices for each neuron
non_neuronal_types_indices <- lapply(non_neuronal_types, function(x){which(sc_ob$neuron_level==x)})

# Sample down to 1k to avoid issues with different numbers of starting cells influencing the overall CV
set.seed(42)
non_neuronal_types_indices <- lapply(non_neuronal_types_indices, sample, 1000)

# Filter to non-neuronal cells
non_neuronal_cells <- subset(sc_ob, cells = unlist(non_neuronal_types_indices))

# Aggregate by tissue type (sum counts within each tissue)
non_neuronal_matrix <- AggregateExpression(
  non_neuronal_cells, 
  group.by = "Tissue",
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE)$RNA

print(paste("Non-neuronal reference matrix dimensions:",
            nrow(non_neuronal_matrix), "x", ncol(non_neuronal_matrix)))
rm(non_neuronal_cells)

# Prepare neuron data for on-the-fly subsampling
neuron_types <- c('CAN', 'DB', 'VB', 'AVJ', 'RIG', 'RMD', 'SIA', 'SMB', 'SMD', 'VD_DD')
names(neuron_types) <- neuron_types

# Subset to just neurons
sc_ob_neuron <- subset(sc_ob, subset = neuron_level %in% neuron_types)

subsample_sizes <- c(10, 50, 100, 500)
n_iterations <- 500

## get indices for each neuron type
neuron_indices <- lapply(neuron_types, function(x){which(sc_ob_neuron$neuron_level==x)})

# Sample down to 1k to avoid issues with different numbers of starting cells influencing the overall CV
set.seed(42)
neuron_indices <- lapply(neuron_indices, sample, 1000)

cat("Neuron indices prepared for on-the-fly subsampling\n")
for(neuron_type in names(neuron_indices)) {
  cat(paste("  ", neuron_type, ":", length(neuron_indices[[neuron_type]]), "cells available\n"))
}

# Load bulk data
bulk <- read.table('Data/231109_bsn12_count.txt', header = 1, row.names = 1)
bulk <- bulk[,6:ncol(bulk)]
colnames(bulk) <- str_split_fixed(colnames(bulk), 'bams.', 2)[,2] |>
  gsub(pattern = '.bam', replace = '', x = _)
bulk <- bulk[,str_split_fixed(colnames(bulk), 'r', 2)[,1] %in% neuron_types]
bulk <- bulk[order(rownames(bulk)),]





deconvolution_results <- list()
correlation_results <- list() 

n_cores <- detectCores() - 1  # Use all cores except one
n_cores <- 4
cat(paste("Setting up Mac parallel processing with", n_cores, "cores\n"))

for(sample_name in colnames(bulk)) {
  cat(paste("\n=== Processing sample:", sample_name, "===\n"))
  
  target_neuron_type <- str_split_fixed(sample_name, 'r', 2)[,1]
  cat(paste("Target neuron type:", target_neuron_type, "\n"))
  
  deconvolution_results[[sample_name]] <- list()
  correlation_results[[sample_name]] <- list()  # NEW
  
  bulk_sample <- bulk[, sample_name, drop = FALSE]
  
  for(subsample_size in subsample_sizes) {
    cat(paste("  Subsample size:", subsample_size, "\n"))
    
    subtracted_results_matrix <- NULL
    
    iteration_results <- pbmclapply(1:n_iterations, function(iteration) {
      
      set.seed(42 + iteration)
      
      # Sample cells from the target neuron type
      target_cells <- sample(neuron_indices[[target_neuron_type]], subsample_size)
      
      # Create subsampled reference for target neuron type
      # Extract raw counts and aggregate manually to avoid AggregateExpression issues
      target_cells_data <- GetAssayData(sc_ob_neuron[, target_cells], assay = "RNA", slot = "counts")
      target_reference_vector <- Matrix::rowSums(target_cells_data)
      
      # Create simple reference matrix with just target neuron type
      neuronal_reference <- data.frame(target_reference_vector)
      colnames(neuronal_reference) <- target_neuron_type
      
      common.genes <- intersect(rownames(bulk_sample), rownames(neuronal_reference))
      
      bulk_use <- bulk_sample[common.genes, , drop = FALSE]
      neuronal_use <- neuronal_reference[common.genes, , drop = FALSE]
      non_neuronal_use <- data.frame(non_neuronal_matrix)[common.genes,]
      
      reference_use <- cbind(neuronal_use, non_neuronal_use)
      reference_use <- cpm(reference_use)
      spm <- apply(log1p(reference_use), 1, max_spm)
      spm <- spm[common.genes]
      
      contaminants <- colnames(non_neuronal_use)
      cell_types_matrix <- matrix(c(target_neuron_type, contaminants), nrow = 1)
      colnames(cell_types_matrix) <- c('target', contaminants)
      rownames(cell_types_matrix) <- sample_name
      
      # Initial Deconvolution Estimates
      r <- reference_use[,cell_types_matrix[1,]]
      A <- nnls(A = as.matrix(log1p(r) * spm),
                b = as.matrix(log1p(bulk_use[[1]]) * spm))$x
      A <- A/sum(A)
      names(A) <- colnames(r)
      initial_estimates <- A
      
      # Final Deconvolution with Contamination Correction
      bulk_subtracted <- subtraction(bulk = bulk_use,
                                     reference = reference_use,
                                     cell_types_matrix = cell_types_matrix,
                                     training_matrix = UNN_train,
                                     specificity_weights = spm,
                                     verbose = FALSE)
      
      subtracted_values <- bulk_subtracted[[1]]
      subtracted_mean <- mean(subtracted_values, na.rm = TRUE)
      subtracted_median <- median(subtracted_values, na.rm = TRUE)
      
      # Final estimates on subtracted data
      A_final <- nnls(A = as.matrix(log1p(r) * spm),
                      b = as.matrix(log1p(bulk_subtracted[[1]]) * spm))$x
      A_final <- A_final/sum(A_final)
      names(A_final) <- colnames(r)
      final_estimates <- A_final
      
      # Train AUROC
      b <- bulk_subtracted[,1]
      names(b) <- rownames(bulk_subtracted)
      train_AUROC <- calc_bulk_auc_one_sample(bulk_vector = log1p(b),
                                              sample_name = sample_name,
                                              training_matrix = UNN_train,
                                              training_genes = rownames(UNN_train),
                                              threshold_list = seq(-5,10,0.01),
                                              sep = 'r')
      
      # Test AUROC
      test_AUROC <- calc_bulk_auc_one_sample(bulk_vector = log1p(b),
                                             sample_name = sample_name,
                                             training_matrix = UNN_test,
                                             training_genes = rownames(UNN_test),
                                             threshold_list = seq(-5,10,0.01),
                                             sep = 'r')
      
      # Return results for this iteration
      list(
        initial_estimates = initial_estimates,
        final_estimates = final_estimates,
        train_AUROC = train_AUROC,
        test_AUROC = test_AUROC,
        subtracted_mean = subtracted_mean,
        subtracted_median = subtracted_median,
        subtracted_values = subtracted_values  # NEW: Return the actual subtracted values
      )
    }, mc.cores = n_cores)
    
    subtracted_list <- lapply(iteration_results, function(x) x$subtracted_values)
    
    gene_names <- names(subtracted_list[[1]])
    subtracted_matrix <- do.call(cbind, subtracted_list)
    rownames(subtracted_matrix) <- gene_names
    colnames(subtracted_matrix) <- paste0("iter_", 1:n_iterations)
    
    # Calculate pairwise correlations
    correlation_analysis <- calculate_pairwise_correlations(subtracted_matrix)
    
    # Calculate gene stability metrics
    gene_stability_analysis <- calculate_gene_stability(subtracted_matrix)
    
    # Store
    correlation_results[[sample_name]][[as.character(subsample_size)]] <- list(
      correlation_metrics = correlation_analysis,
      gene_stability_metrics = gene_stability_analysis
    )
    
    # Clean up clean up *sing this part*
    for(i in 1:length(iteration_results)) {
      iteration_results[[i]]$subtracted_values <- NULL
    }
    
    # Store results for this subsample size
    deconvolution_results[[sample_name]][[as.character(subsample_size)]] <- iteration_results
    
    # Additional cleanup after each subsample size *so much mess to clean...*
    cat(paste("    Completed subsample size", subsample_size, "- cleaning memory\n"))
    rm(subtracted_matrix, subtracted_list, correlation_analysis)  # NEW: Clean up correlation data
    gc(verbose = FALSE)
  }
  
  # Final cleanup after each sample
  cat(paste("Completed sample", sample_name, "- performing garbage collection\n"))
  gc(verbose = FALSE)
}


cv_results <- calculate_cv_results(deconvolution_results, correlation_results)


cv_summary <- create_cv_summary_table(cv_results)
correlation_by_size <- aggregate(mean_pairwise_correlation ~ subsample_size,
                                 data = cv_summary, FUN = mean, na.rm = TRUE)
cv_by_size <- aggregate(cbind(final_purity_cv, test_auroc_cv, subtracted_mean_cv, subtracted_median_cv) ~ subsample_size,
                        data = cv_summary, FUN = mean, na.rm = TRUE)

means_by_size <- aggregate(cbind(final_purity_mean, test_auroc_mean, subtracted_mean_mean, subtracted_median_mean) ~ subsample_size,
                           data = cv_summary, FUN = mean, na.rm = TRUE)
print(means_by_size)

# Save results
saveRDS(list(
  deconvolution_results = deconvolution_results,
  correlation_results = correlation_results, 
  cv_results = cv_results,
  cv_summary = cv_summary
), "sample_wise_deconvolution_results_with_correlations.rds")


sample_results <- readRDS('sample_wise_deconvolution_results_with_correlations.rds')

cv_summary <- sample_results$cv_summary

## Make plots




colnames(cv_summary)
cv_summary |> ggplot() +
  geom_boxplot(aes(x = as.factor(subsample_size),
                   y = test_auroc_mean,
                   fill = as.factor(subsample_size),),
               position = 'dodge2', notch = T)

cv_summary |> ggplot() +
  geom_boxplot(aes(x = log10(subsample_size),
                   y = mean_pairwise_correlation,
                   fill = as.factor(subsample_size),),
               position = 'dodge2', notch = T)

cv_summary |> ggplot() +
  geom_boxplot(aes(x = log10(subsample_size),
                   y = final_purity_mean,
                   fill = as.factor(subsample_size),),
               position = 'dodge2', notch = T)

cv_summary |> ggplot() +
  geom_boxplot(aes(x = log10(subsample_size),
                   y = final_purity_cv,
                   fill = as.factor(subsample_size),),
               position = 'dodge2', notch = T)

cv_summary |> ggplot() +
  geom_boxplot(aes(x = log10(subsample_size),
                   y = subtracted_median_mean,
                   fill = as.factor(subsample_size),),
               position = 'dodge2', notch = T)


summary(glm(mean_gene_cv ~ log10(subsample_size), data = cv_summary))
summary(glm(final_purity_mean ~ log10(subsample_size), data = cv_summary))
summary(glm(test_auroc_mean ~ log10(subsample_size), data = cv_summary))
summary(glm(mean_pairwise_correlation ~ log10(subsample_size), data = cv_summary))


results <- cv_summary |>
  group_by(neuron_type) |>
  group_modify(~{
    model <- glmmTMB(log10(subtracted_mean_mean) ~ log10(subsample_size) + (1|sample),
                     family = gaussian(link='identity'),
                     data = .x)
    
    tibble(
      slope = fixef(model)$cond["log10(subsample_size)"],
      slope_se = sqrt(vcov(model)$cond["log10(subsample_size)", "log10(subsample_size)"]),
      t_value = slope / slope_se,
      p_value = summary(model)$coefficients$cond["log10(subsample_size)", "Pr(>|z|)"]
    )
  })

results |> select(neuron_type, slope, slope_se, t_value, p_value)



cv_summary |> 
  ggplot(aes(x = subsample_size,
             y = initial_purity_mean),) +
  geom_violin(aes(fill = as.factor(subsample_size)), ) +
  geom_smooth(
    method = 'lm', se = F, color = 'black') +
  scale_x_continuous(transform = 'log10') +
  ylab('Average neuron composition estimate per unaltered sample') +
  xlab('single cell reference cluster size') +
  ggtitle('LittleBites reference depth\nvs\nneuron composition estimate') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = 'black'))
ggsave('figures/Figure_2/reference_subsampling_initial_composition.pdf', width = 10, height = 10)

cv_summary |> 
  ggplot(aes(x = subsample_size,
             y = initisubtracted_mean_mean),) +
  geom_violin(aes(fill = as.factor(subsample_size)), ) +
  geom_smooth(aes(group = neuron_type, color = neuron_type),
    method = 'lm', se = F) +
  scale_x_continuous(transform = 'log10') +
  scale_y_continuous(transform = 'log10') +
  ylab('Cleaned sample mean expression') +
  xlab('single cell reference cluster size') +
  ggtitle('LittleBites reference depth\nvs\nsubtracted_sample_') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = 'black'))
ggsave('figures/Figure_2/reference_subsampling_initial_composition.pdf', width = 10, height = 10)

cv_summary |> 
  ggplot(aes(x = subsample_size,
             y = mean_gene_cv),) +
  geom_violin(aes(fill = as.factor(subsample_size)), ) +
  geom_smooth(
    method = 'lm', se = F) +
  scale_x_continuous(transform = 'log10') +
  ylab('Average Coefficient of Variation per gene') +
  xlab('single cell reference cluster size') +
  ggtitle('LittleBites reference depth\nvs\ngene expression stability') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = 'black'))
ggsave('figures/Figure_2/reference_subsampling_gene_CV.pdf', width = 10, height = 10)

cv_summary |> 
  ggplot(aes(x = subsample_size,
             y = final_purity_mean),) +
  geom_violin(aes(fill = as.factor(subsample_size)), ) +
  geom_smooth(
    method = 'lm', se = F) +
  scale_x_continuous(transform = 'log10') +
  ylab('Average sample purity estimate (NNLS model) after subtraction') +
  xlab('single cell reference cluster size') +
  ggtitle('LittleBites reference depth\nvs\npost-subtraction purity estimate') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = 'black'))
ggsave('figures/Figure_2/reference_subsampling_subtracted_purity.pdf', width = 10, height = 10)

cv_summary |> 
  ggplot(aes(x = subsample_size,
             y = test_auroc_mean),) +
  geom_violin(aes(fill = as.factor(subsample_size)), ) +
  geom_smooth(
    method = 'lm', se = F) +
  scale_x_continuous(transform = 'log10') +
  ylab('Average AUROC per sample') +
  xlab('single cell reference cluster size') +
  ggtitle('LittleBites reference depth\nvs\nground-truth AUROC') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = 'black'))
ggsave('figures/Figure_2/reference_subsampling_AUROC.pdf', width = 10, height = 10)

cv_summary |> 
  ggplot(aes(x = subsample_size,
             y = subtracted_mean_cv),) +
  geom_violin(aes(fill = as.factor(subsample_size)), ) +
  geom_smooth(
    method = 'lm', se = F) +
  scale_x_continuous(transform = 'log10') +
  ylab('CV of cleaned sample mean expression') +
  xlab('single cell reference cluster size') +
  ggtitle('LittleBites reference depth\nvs\nsubtracted mean value variance') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = 'black'))
ggsave('figures/Figure_2/reference_subsampling_subtracted_mean_CV.pdf', width = 10, height = 10)


