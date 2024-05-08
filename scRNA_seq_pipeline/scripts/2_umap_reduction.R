library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)

principal_components <- snakemake@params[[8]]

# load in the pruned data after performing QC
Sample1 = readRDS(file=snakemake@params[[1]])
Sample2 = readRDS(file=snakemake@params[[2]])

name_of_sample_1 <- snakemake@params[[6]]
name_of_sample_2 <- snakemake@params[[7]]

#normalize data and find variable features
#Sample1$Sample2<-'Sample1'
Sample1$Sample2 <- name_of_sample_1
Sample1<-NormalizeData(Sample1, verbose=FALSE)
Sample1<-FindVariableFeatures(Sample1,selection.method = "vst", nfeatures = 2000)

#Sample2$Sample2<-'Sample2'
Sample2$Sample2 <- name_of_sample_2
Sample2<-NormalizeData(Sample2, verbose=FALSE)
Sample2<-FindVariableFeatures(Sample2, selection.method="vst", nfeatures=2000)

#FindIntegrationAnchors function takes a list of Seurat objects as input and identifies anchors and then uses these anchors to integrate two datasets together with integrate data, trying here to use 2 datasets
STRESS.anchors<-FindIntegrationAnchors(object.list = list(Sample1,Sample2), dims= 1:20)
STRESS.combined<-IntegrateData(anchorset=STRESS.anchors, dims=1:20)

#Create a violin plot of features of interest
pdf(file=snakemake@output[[1]])
VlnPlot(STRESS.combined, features=c(snakemake@params[[3]]), split.by="Sample2", assay="RNA", combine=FALSE)
#VlnPlot(STRESS.combined, features=c(snakemake@params[[3]]), split.by=name_of_sample_2, assay="RNA", combine=FALSE)
dev.off()

#Perform integrated analysis, here can run a single integrated analysis on all cells 
DefaultAssay(STRESS.combined)<-"integrated"

#Run the standard workflow for visualization and clustering 
STRESS.combined<-ScaleData(STRESS.combined, verbose=FALSE)
STRESS.combined<-RunPCA(STRESS.combined, npcs = 30, verbose = FALSE)

####User PC caclulation####
pct <- STRESS.combined[["pca"]]@stdev / sum(STRESS.combined[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
# last point where change of % of variation is more than 0.1%.
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

pcs <- min(co1,co2)

if (principal_components != 'automatic'){
pcs <- strtoi(principal_components)
}

p1 <- ElbowPlot(STRESS.combined, ndims = 30)

elbow_plot_title <- paste("Calculated number of pcs:", as.character(pcs))
pdf(file=snakemake@output[[5]])
p1 + labs(title= elbow_plot_title)+ theme(plot.title = element_text(hjust=0.5))
dev.off()

#t-SNE and Clustering
####USER INPUT 2####
STRESS.combined <- RunUMAP(STRESS.combined, reduction = "pca", dims = 1:pcs)
STRESS.combined <- FindNeighbors(STRESS.combined, reduction = "pca", dims = 1:pcs)
STRESS.combined <- FindClusters(STRESS.combined, resolution = snakemake@params[[5]])

# umap Visualization
pdf(file=snakemake@output[[2]])
p2 <- DimPlot(STRESS.combined, reduction = "umap", group.by = "Sample2")
#p2 <- DimPlot(STRESS.combined, reduction = "umap", group.by = name_of_sample_2)
p3 <- DimPlot(STRESS.combined, reduction = "umap", label = TRUE)
plot_grid(p3)

dev.off()

pdf(file=snakemake@output[[4]], width = 16)
VlnPlot(STRESS.combined, features = c(snakemake@params[[4]]) , split.by="Sample2", assay="RNA",  combine=TRUE, cols = c("grey", "red"))
#VlnPlot(STRESS.combined, features = c(snakemake@params[[4]]) , split.by=name_of_sample_2, assay="RNA",  combine=TRUE, cols = c("grey", "red"))
dev.off()

saveRDS(STRESS.combined, file = snakemake@output[[3]])