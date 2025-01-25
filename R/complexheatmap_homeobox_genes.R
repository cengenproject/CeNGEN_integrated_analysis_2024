library(pheatmap)
library(dplyr)
library(stringr)
library(grid)
library(ComplexHeatmap)
library(wbData)

## pick wormbase reference

ws289 <- wbData::wb_load_gene_ids('289')


### load in data
integrated <- read.csv("Data/Integrated_thresholded/240721_Integrated_bsn12_cpm_unthresholded.csv.gz")
subtracted <- read.csv('Data/bsn12_bulk_subtracted_TMM_051624.tsv.gz', sep = '\t')
bulk <- read.csv('Data/bsn12_bulk_TMM_051624.tsv.gz', sep = '\t')
prop2count <- read.csv('Data/sc_biorep_prop2count_050124.tsv.gz', sep = '\t')


## aggregate replicates
aggr_prop2count <- prop2count
colnames(aggr_prop2count) <-str_split_fixed(colnames(aggr_prop2count),"__",2)[,1]
aggr_prop2count <- data.frame(vapply(unique(colnames(aggr_prop2count)), function(x)
  rowMeans(aggr_prop2count[,colnames(aggr_prop2count)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_prop2count)) ))
aggr_prop2count$VD <- aggr_prop2count$VD_DD
aggr_prop2count$DD <- aggr_prop2count$VD_DD

aggr_bulk <- bulk
colnames(aggr_bulk) <-str_split_fixed(colnames(aggr_bulk),"r",2)[,1]
aggr_bulk <- data.frame(vapply(unique(colnames(aggr_bulk)), function(x)
  rowMeans(aggr_bulk[,colnames(aggr_bulk)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_bulk)) ))

aggr_subtracted <- subtracted
colnames(aggr_subtracted) <-str_split_fixed(colnames(aggr_subtracted),"r",2)[,1]
aggr_subtracted <- data.frame(vapply(unique(colnames(aggr_subtracted)), function(x)
  rowMeans(aggr_subtracted[,colnames(aggr_subtracted)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_subtracted)) ))


## subset to common neurons and genes
neurons <- table(c(colnames(integrated), colnames(aggr_subtracted), colnames(aggr_prop2count))) |> 
  {\(x) x[x==max(x)]}() |>
  names()

genes <- table(c(rownames(integrated), rownames(aggr_subtracted), rownames(aggr_prop2count))) |> 
  {\(x) x[x==max(x)]}() |>
  names()

integrated <- integrated[genes, neurons] |> log1p()
aggr_bulk <- aggr_bulk[genes, neurons] |> log1p()
aggr_subtracted <- aggr_subtracted[genes, neurons] |> log1p()
aggr_prop2count <- aggr_prop2count[genes, neurons] |> log1p()



## load in homeobox gene protein level annotations
homeobox_genes <- readRDS('references/072820_homeobox_truth.rds')

## collapse homeobox gene expression to bulk cell annotations
homeobox_genes$AWC <- homeobox_genes$AWC_OFF + homeobox_genes$AWC_ON
homeobox_genes$IL2 <- homeobox_genes$IL2_DV + homeobox_genes$IL2_LR
homeobox_genes$DA <- homeobox_genes$DA + homeobox_genes$DA9
homeobox_genes$DB <- homeobox_genes$DB + homeobox_genes$DB01
homeobox_genes$RMD <- homeobox_genes$RMD_DV + homeobox_genes$RMD_LR
homeobox_genes$RME <- homeobox_genes$RME_DV + homeobox_genes$RME_LR
homeobox_genes$VA <- homeobox_genes$VA + homeobox_genes$VA12
homeobox_genes$VB <- homeobox_genes$VB + homeobox_genes$VB01 + homeobox_genes$VB02

homeobox_genes[homeobox_genes>0] <- 1

## subset to common neurons
homeobox_genes <- homeobox_genes[,neurons]

## just include homeobox genes
integrated <- integrated[rownames(homeobox_genes),]
aggr_bulk <- aggr_bulk[rownames(homeobox_genes),]
aggr_subtracted <- aggr_subtracted[rownames(homeobox_genes),]
aggr_prop2count <- aggr_prop2count[rownames(homeobox_genes),]


rownames(integrated) <- i2s(rownames(integrated), ws289)
rownames(aggr_bulk) <- i2s(rownames(aggr_bulk), ws289)
rownames(aggr_subtracted) <- i2s(rownames(aggr_subtracted), ws289)
rownames(aggr_prop2count) <- i2s(rownames(aggr_prop2count), ws289)

rownames(homeobox_genes) <- i2s(rownames(homeobox_genes), ws289)


pheatmap(integrated, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("white", "orange", "maroon", "navy"))(10),
         fontsize_row = 5, fontsize_col = 5, cellheight = 5, cellwidth = 5, 
         breaks = c(0:10))


draw_overlay <- function(j, i, x, y, width, height, fill) {
  value <- homeobox_genes[i, j]  # Get the binary value from df_2
  square_color <- ifelse(value == 1, "black", "white")  # Define color
  grid.rect(x = x, y = y, width = unit(2, "mm"), height = unit(2, "mm"),
            gp = gpar(fill = square_color, col = NA))
}

p2c_hm <- Heatmap(aggr_prop2count,
  name = " ",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = colorRampPalette(c("white", "orange", "maroon", "navy"))(20),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA))
    draw_overlay(j, i, x, y, width, height, fill)
  }
)

int_hm <- Heatmap(integrated,
                  name = " ",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  col = colorRampPalette(c("white", "orange", "maroon", "navy"))(20),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA))
                    draw_overlay(j, i, x, y, width, height, fill)
                  }
)

bulk_hm <- Heatmap(aggr_bulk,
                  name = " ",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  col = colorRampPalette(c("white", "orange", "maroon", "navy"))(20),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA))
                    draw_overlay(j, i, x, y, width, height, fill)
                  }
)

sub_hm <- Heatmap(aggr_subtracted,
                  name = " ",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  col = colorRampPalette(c("white", "orange", "maroon", "navy"))(20),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA))
                    draw_overlay(j, i, x, y, width, height, fill)
                  }
)


n_rows <- nrow(homeobox_genes)
n_cols <- ncol(homeobox_genes)
cell_size <- .2  

heatmap_width <- cell_size * n_cols
heatmap_height <- cell_size * n_rows

pdf("prop2count_homeobox_GT_250124.pdf", width = heatmap_width, height = heatmap_height)
draw(p2c_hm, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
pdf("integrated_homeobox_GT_250124.pdf", width = heatmap_width, height = heatmap_height)
draw(int_hm, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
pdf("bulk_homeobox_GT_250124.pdf", width = heatmap_width, height = heatmap_height)
draw(bulk_hm, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
pdf("subtracted_homeobox_GT_250124.pdf", width = heatmap_width, height = heatmap_height)
draw(sub_hm, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
