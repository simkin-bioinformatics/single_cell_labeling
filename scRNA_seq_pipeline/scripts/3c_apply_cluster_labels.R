#rename clusters. Use the marker genes known to associate with cell types for this
#the number of new names needs to match the number of clusters, and duplicates are allowed
#example given below with new.cluster.ids

library(Seurat)
library(dplyr)
library(Matrix)
library(reticulate)
library(ggplot2)
# install.packages('BiocManager')
# BiocManager::install('limma')

#setwd('D:/yhtze/Merfish')

options(scipen = 100)

STRESS.combined = readRDS(file=snakemake@input[[2]])

Idents(STRESS.combined) <- STRESS.combined@meta.data$seurat_clusters

STRESS.combined@meta.data$seurat_clusters

cluster_labels <- read.table(file = snakemake@input[[1]],sep="\t", header = TRUE, row.names=1)

new.cluster.ids <- cluster_labels$cluster_names
length(new.cluster.ids)
levels(STRESS.combined)
names(new.cluster.ids) <- levels(STRESS.combined)
STRESS.combined <- RenameIdents(STRESS.combined, new.cluster.ids)

pdf(file=snakemake@output[[2]], width=24, height=20)
DimPlot(STRESS.combined, reduction = "umap", label=TRUE)
dev.off()

saveRDS(STRESS.combined,file=snakemake@output[[1]])
