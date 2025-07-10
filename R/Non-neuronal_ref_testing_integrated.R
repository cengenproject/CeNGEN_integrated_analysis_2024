#### non-neuronal testing in integrated at single cell data

library(ggplot2)



non_neuronal_genes <- read.table('references/ubituiqtous_and_nonNeuronal_gt_genes_matrix_042222.tsv')

non_neuronal_genes <- non_neuronal_genes[rowSums(non_neuronal_genes) == 0, ] |> rownames()
non_neuronal_genes



cengen_sc_2_bulk
int_2 <- read.table('Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_2_063025.csv', sep = ',')
int_2 <- (int_2 > 0)*1

int_4 <- read.table('Data_out/Integrated_thresholded/Integrated_bsn12_cpm_threshold_4_063025.csv', sep = ',')
int_4 <- (int_4 > 0)*1


cengen_sc_2_bulk
neurons <- intersect(colnames(cengen_sc_2_bulk), colnames(int_4))

non_neuronal_genes_check <- intersect(non_neuronal_genes, rownames(int_2)) |> 
  intersect(rownames(cengen_sc_2_bulk)) |>
  sort()

blank <- rep(0, length(intersect(non_neuronal_genes, rownames(int_2))))
names(blank) <- intersect(non_neuronal_genes, rownames(int_2))




sum_detected_ <- rowSums(int_2[non_neuronal_genes_check,neurons])
sum_detected4_ <- rowSums(int_4[non_neuronal_genes_check,neurons])
sum_detected_sc2_ <- rowSums(cengen_sc_2_bulk[non_neuronal_genes_check,neurons] > 0)
sum_detected_sc4_ <- rowSums(cengen_sc_4_bulk[non_neuronal_genes_check,neurons] > 0)

sum_detected <- blank
sum_detected4 <- blank
sum_detected_sc2 <-blank
sum_detected_sc4 <- blank

sum_detected[non_neuronal_genes_check] <- sum_detected_
sum_detected4[non_neuronal_genes_check] <- sum_detected4_
sum_detected_sc2[non_neuronal_genes_check] <- sum_detected_sc2_
sum_detected_sc4[non_neuronal_genes_check] <- sum_detected_sc4_

dim(cengen_sc_4_bulk)

hist(sum_detected4)

summary(sum_detected_sc4)
summary(sum_detected4)
summary(sum_detected_sc2)
summary(sum_detected)

collected <- data.frame(gene = rep(names(sum_detected), 4),
           value = c(sum_detected,
                     sum_detected_sc2,
                     sum_detected4,
                     sum_detected_sc4),
           source = c(rep('int', length(sum_detected)),
                      rep('sc', length(sum_detected)),
                      rep('int', length(sum_detected)),
                      rep('sc', length(sum_detected))),
           FDR_threshold = c(rep('14 percent', length(sum_detected)*2),
                         rep('8.9 percent', length(sum_detected)*2)))

ggplot(collected) + 
  geom_violin(aes(x = threshold, y = value, fill = source),
               position = 'dodge', scale = 'area')
ggplot(collected) + 
  geom_boxplot(aes(x = FDR_threshold, y = value, fill = source),
              position = 'dodge', notch = T) +
  theme_classic(base_size = 10) +
  ylab('ncells detected') +
  xlab('FDR threshold') +
  ggtitle('Detection rate, per cell\nnon-neuronal markers') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave('figures/non_neuronal_gene_detection_thresolds.pdf')
plot(sum_detected_sc4, sum_detected4)
abline(0,1)
sum_detected[sum_detected>50]
sum(sum_detected4 > 10)



