####QC####
setwd("~/Desktop/Data/scRNA_Seq/Tapper_Stress_Fall_2022/")
#To do
#1 Annotate Cell Types, dim=20
#2 Subset Neurons
#3 Look for cFos expression across timepoints
#4 Look for markers in cFos+ cluster across all neurons
library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)


T0<- Read10X(data.dir = "T0/filtered_feature_bc_matrix/" )
T0<- CreateSeuratObject(counts = T0, project = "CTL")
T0


T5<- Read10X(data.dir = "T5/filtered_feature_bc_matrix/")
T5<- CreateSeuratObject(counts=T5, project= "STRESS")
T5

T0$T5<-"T0"

####USER INPUT 1####
T0[["percent.mt"]] <- PercentageFeatureSet(T0, pattern = "^mt-")
VlnPlot(T0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(T0, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(T0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

T0<-subset(T0, subset=nFeature_RNA>200 & nFeature_RNA < 7000 & percent.mt < 10)
T0<-NormalizeData(T0, verbose=FALSE)
T0<-FindVariableFeatures(T0,selection.method = "vst", nfeatures = 2000)

T5$T5<-"T5"
T5[["percent.mt"]] <- PercentageFeatureSet(T5, pattern = "^mt-")
VlnPlot(T5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(T5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(T5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

T5<-subset(T5, subset=nFeature_RNA>200 & nFeature_RNA < 7000 & percent.mt < 10)
T5<-NormalizeData(T5, verbose=FALSE)
T5<-FindVariableFeatures(T5, selection.method="vst", nfeatures=2000)


#FindIntegrationAnchors function takes a list of Seurat objects as input and identifies anchors and then uses these anchors to integrate two datasets together with integrate data, trying here to use 2 datasets
STRESS.anchors<-FindIntegrationAnchors(object.list = list( T0,T5), dims= 1:20)
STRESS.combined<-IntegrateData(anchorset=STRESS.anchors, dims=1:20)

VlnPlot(STRESS.combined, features=c("Fos", "Map2", "Gad2"), split.by="T5", assay="RNA", combine=FALSE)
#Perform integrated analysis, here can run a single integrated analysis on all cells 
DefaultAssay(STRESS.combined)<-"integrated"

#Run the standard workflow for viauslization and clustering 
STRESS.combined<-ScaleData(STRESS.combined, verbose=FALSE)
STRESS.combined<-RunPCA(STRESS.combined, npcs = 30, verbose = FALSE)


####User PC caclulation####
pct <- STRESS.combined[["pca"]]@stdev / sum(STRESS.combined[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
#t-SNE and Clustering
####USER INPUT 2####
STRESS.combined <- RunUMAP(STRESS.combined, reduction = "pca", dims = 1:24)
STRESS.combined <- FindNeighbors(STRESS.combined, reduction = "pca", dims = 1:24)
STRESS.combined <- FindClusters(STRESS.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(STRESS.combined, reduction = "umap", group.by = "T5")
p2 <- DimPlot(STRESS.combined, reduction = "umap", label = TRUE)
plot_grid(p2)


DimPlot(STRESS.combined, reduction = "umap", label=TRUE)

VlnPlot(STRESS.combined, features = c("Map2", "Thy1","Sst", "Gad1", "Gad2", "Fos") , split.by="T5", assay="RNA",  combine=TRUE,
        cols = c("grey", "red"))

####USER INPUT 3####
####Annotation####

cluster.markers <- FindAllMarkers(STRESS.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(cluster.markers, file="IPN_Stress_Cluster_Markers.csv")
write.table(cluster.markers, file='T5_vs_T0_cluster_markers.tsv', quote=FALSE, sep='\t', col.names = NA)

STRESS.combined <- RenameIdents(STRESS.combined,  `0` = "Neuron", `1` = "Oligo", `2` = "Microglia", 
                                `3` = "Astrocyte", `4` = "Oligo", `5` = "Endothelial", `6` = "OPC", `7` = "Microglia", `8` = "Oligodendrocyte", `9` = "Neuron", 
                                `10` = "Neuron", `11` = "Oligo", `12` = "Endothelial", `13` = "Mural", `14` = "14", `15` = "Neuron",
                                `16` = "16")


####Subset Neurons####


STRESS.combined_C2<-subset(STRESS.combined, idents = c("0", "6","10","12", "18",))
STRESS.combined_C2<-subset(STRESS.combined, idents = c("Neuron"))


STRESS.combined_C2<-ScaleData(STRESS.combined_C2)
STRESS.combined_C2<-NormalizeData(STRESS.combined_C2)
FeaturePlot(STRESS.combined_C2, "Sst")+(ggplot2::scale_fill_continuous(limits = c(0.0,1.0),, breaks = c(0.0, 0.5, 1.0)))

pct <- STRESS.combined[["pca"]]@stdev / sum(STRESS.combined_C2[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2







ElbowPlot(STRESS.combined_C2)
STRESS.combined_C2 <- FindNeighbors(STRESS.combined_C2, graph.name = "Test", dims = 1:20)
STRESS.combined_C2 <- FindClusters(STRESS.combined_C2, graph.name = "Test", resolution = 0.5, algorithm = 1, verbose = TRUE)
STRESS.combined_C2 <- RunUMAP(STRESS.combined_C2, dims = 1:20)


DefaultAssay(STRESS.combined_C2)<-"RNA"
VlnPlot(STRESS.combined_C2,features=c( "Sst","Fos","Gad2", "Gad1","Htr2c", "Grm1"), assay="RNA", pt.size=0, split.by = "T5")


####USER INPUT 4####
####Differential Expression####

STRESS.combined$celltype.T5<-paste(Idents(STRESS.combined), STRESS.combined$T5, sep="_")
STRESS.combined$celltype<-Idents(STRESS.combined)
Idents(STRESS.combined)<-"celltype.T5"


Neuron_stress<-FindMarkers(STRESS.combined, ident.1 = "19_T5", ident.2 = "19_T0", verbose=FALSE)
write.csv(Neuron_stress, file = "19_T5_vs_T0.csv")

neuron_data<-read.csv("18_T5_vs_T0.csv", stringsAsFactors = TRUE)  
names(neuron_data)[1]<-as.character('gene')
neuron_plot<-ggplot(data=neuron_data, aes(x=avg_log2FC, y=-log10(p_val)))+geom_point()+theme_minimal()
p2<-neuron_plot+geom_vline(xintercept = c(-0.6, 0.6, col="red")+geom_hline(yintercept=-log10(0.05),col="red"))
neuron_data$diffexpressed="NO"
neuron_data$diffexpressed[neuron_data$avg_log2FC > 0.5849 & neuron_data$p_val< 0.05] <- "UP"
neuron_data$diffexpressed[neuron_data$avg_log2FC < -0.5849 & neuron_data$p_val< 0.05] <- "DOWN"
neuron_plot<-ggplot(data=neuron_data, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed)) + geom_point() + theme_minimal()
p2<- neuron_plot + geom_vline(xintercept = c(-0.6, 0.6), col="red")+ geom_hline(yintercept= -log10(0.05), col="red")
mycolors<-c("blue", "red", "black")
names(mycolors)<-c("DOWN", "UP", "NO")
p3<-p2 +scale_color_manual(values=mycolors)
neuron_data$delabel<-NA
neuron_data$delabel[neuron_data$diffexpressed != "NO"]<-as.character(neuron_data$gene)[neuron_data$diffexpressed != "NO"]
p4<-ggplot(data=neuron_data, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel))+
  geom_point()+theme_minimal()+geom_text()

write.csv(neuron_data, file = "18_T5_vs_T0.csv")

#neuron_data<-read.csv("15_T5_vs_T0.csv", stringsAsFactors = TRUE)

library(ggrepel)
p5<-ggplot(data= neuron_data, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel))+
  geom_point()+
  theme_minimal()+
  ggrepel::geom_text_repel()+
  scale_color_manual(values=c("blue", "black", "red"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="red")+
  geom_hline(yintercept = -log10(0.05), col="red")
p5 + labs(title= "Type II Spiral Ganglion T5 vs T0")+ theme(plot.title = element_text(hjust=0.5))
print(p5)
write.csv(neuron_data, file = "18_T5_vs_T0.csv")



####Counts per Cell####

DefaultAssay(STRESS.combined_C2)<-"RNA"
new5_t5<-subset(STRESS.combined_C2, idents = c("5_T5"))


write.table(as.matrix(GetAssayData(object = new5_t5, slot = "data")), 
            'new5_T5.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)




