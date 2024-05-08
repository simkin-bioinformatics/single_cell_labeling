####USER INPUT 3####
####Annotation####
library(dplyr)
library(Seurat)
library(cowplot)

STRESS.combined = readRDS(file=snakemake@input[[1]])


cluster.markers <- FindAllMarkers(STRESS.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.table(cluster.markers, file=snakemake@output[[1]], quote=FALSE, sep='\t', col.names = NA)