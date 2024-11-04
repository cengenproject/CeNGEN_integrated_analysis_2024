## libraries ----
library(pbapply)
library(dplyr)
library(wbData)
library(tidyverse)
library(reshape)
library(ggplot2)
library(ggrepel)
library(ggrastr)

## load the data ----


sc_biorep <- readRDS('Data/CeNGEN_biorep_aggregate_counts_data_240709.RDS')


biorep_counts_melt <- melt(sc_biorep$counts)

biorep_proportions_melt <- melt(sc_biorep$proportions)

biorep_cell_number <- data.frame(row.names = names(sc_biorep$nCells),
                nCell = sc_biorep$nCells)


melted_merged <- data.frame(gene = biorep_counts_melt$X1,
                            replicate = biorep_counts_melt$X2,
                            counts = biorep_counts_melt$value,
                            proportions = biorep_proportions_melt$value,
                            cell_number = biorep_cell_number[biorep_counts_melt$X2, 'nCell.Freq'])


## log transform counts, and remove 0 and 1 proportions to eliminate infinite values then calculate the logit ----

melted_merged_log <- melted_merged[melted_merged$proportions > 0 & melted_merged$proportions < 1,]
melted_merged_log$log_counts <- log(melted_merged_log$counts)
melted_merged_log$log_cell_number <- log(melted_merged_log$cell_number)

melted_merged_log$logit_proportions <- log((melted_merged_log$proportions/(1-melted_merged_log$proportions)))

melted_merged_log$odds_ratio <- melted_merged_log$proportions/(1-melted_merged_log$proportions)


repl <- 'SIA__unc-47_2'

melted_merged_log  |> filter(replicate == repl) |>
  arrange(proportions) |>
  ggplot() +
  geom_point_rast(aes(x = log(odds_ratio) + log(cell_number), y = log(counts)),
                  alpha = 0.5, raster.dpi = 300) +
  geom_abline(slope = 1, intercept = 0, color='red') +
  labs(x = 'log(P/(1-P)) + log(nCell)', y = 'log(counts)') +
  theme_classic(base_size =  20) +
  theme(legend.position = '', 
        axis.text = element_text(color = 'black', face = 'bold'),
        axis.title = element_text(color = 'black', face = 'bold'))
ggsave('figures/single_replicate_nbModel_real_counts_SIA__unc-47_2_240531.pdf', width = 5, height = 5)


## define function for simulation

sim_1_nb <- function(N=300,success = 1, prob = .9){
  y <- rnbinom(N, success, prob)
  
  data.frame(C = sum(y),
             P = mean(y > 0))
}


# set.seed(123)
# res_nb_var <- tibble(N = rpois(2000, 10000),
#                      s = rpois(2000, 10)/10,
#                      p = rbeta(2000, 0.5, 0.5)) |>
#   mutate(res = pmap(list(N,s,p), sim_1_nb)) |>
#   unnest(res)



## simulate a cluster with 1000 cells and 2000 genes
set.seed(123)
res_nb_var <- tibble(N = 1000,
                     s = rpois(1000, 10)/10,
                     p = rbeta(1000, 0.5, 0.5)) |>
  mutate(res = pmap(list(N,s,p), sim_1_nb)) |>
  unnest(res)

res_nb_var <- res_nb_var[res_nb_var$P < 1 & res_nb_var$P > 0,]

ggplot(res_nb_var) +
  geom_abline(slope = 1, intercept = 0, color='red') +
  geom_point_rast(aes(x = log(P/(1-P))+log(N), y=log(C)),
                  alpha = 0.5, raster.dpi = 300) +
  labs(x = 'log(P/(1-P)) + log(nCell)', y = 'log(counts)') +
  theme_classic(base_size =  20) +
  theme(legend.position = '', 
        axis.text = element_text(color = 'black', face = 'bold'),
        axis.title = element_text(color = 'black', face = 'bold')) 
ggsave('figures/single_replicate_nbModel_simulated_counts_n10000_240531.pdf', width = 5, height = 5)

