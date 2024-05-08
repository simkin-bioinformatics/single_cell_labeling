library(dplyr)
library(Seurat)
library(cowplot)

# Load the dataset
sample.data <- Read10X(data.dir = snakemake@input[[1]])

# Initialize the Seurat object with the raw (non-normalized data).
sample <- CreateSeuratObject(counts = sample.data, project = snakemake@params[[1]])
sample

# add mitochondrial count column
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt-")

# Visualize QC metrics
pdf(file=snakemake@output[[1]])
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# eliminate parts of the sample based on number of features and percent mt
sample<-subset(sample, subset=nFeature_RNA>snakemake@params[[2]] & nFeature_RNA < snakemake@params[[3]] & percent.mt < snakemake@params[[4]])

saveRDS(sample, file = snakemake@output[[2]])