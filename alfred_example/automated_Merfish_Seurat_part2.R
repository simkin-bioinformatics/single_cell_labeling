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

OB_s=readRDS(file = "OB_s_object.rds")

Idents(OB_s) <- OB_s@meta.data$seurat_clusters

OB_s@meta.data$seurat_clusters

cluster_labels <- read.table(file = 'actual_cluster_labels.txt',sep="\t", header = TRUE, row.names=1)

new.cluster.ids <- cluster_labels$cluster_names
length(new.cluster.ids)
levels(OB_s)
names(new.cluster.ids) <- levels(OB_s)
OB_s <- RenameIdents(OB_s, new.cluster.ids)

DimPlot(OB_s, reduction = "umap", label = TRUE, label.size=6)#+theme(legend.key.size = unit(0.01, 'cm'))+theme(legend.position="bottom", legend.text = element_text(size=5))
#output data file
write.csv(OB_s@active.ident,'unknown_thing.csv')

######
# Extra Plots I've written

###plot xy of each cluster
dfall<- data.frame(OB_s$center.x[WhichCells(object = subset(OB_s,subset = sample =="OB1"))],OB_s$center.y[WhichCells(object = subset(OB_s,subset = sample =="OB1"))],OB_s@active.ident[WhichCells(object = subset(OB_s,subset = sample =="OB1"))])
x <- OB_s$center.x[WhichCells(object = subset(OB_s,subset = sample =="OB1"))]
y <- OB_s$center.y[WhichCells(object = subset(OB_s,subset = sample =="OB1"))]
group <- OB_s@active.ident[WhichCells(object = subset(OB_s,subset = sample =="OB1"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=0.2)+ theme_linedraw()+coord_fixed(ratio = 1)+theme(legend.key.size = unit(0.01, 'cm'))+theme(legend.position="bottom", legend.text = element_text(size=5))
ggsave('file1.pdf')
###plot xy of each cluster in subregion
xmin <- 5000
xmax <- 6250
ymin <- 3750
ymax <- 5000
dfall<- data.frame(OB_s$center.x[WhichCells(object = subset(OB_s,subset = sample =="OB1" & center.x < xmax & center.x > xmin & center.y > ymin & center.y < ymax))],OB_s$center.y[WhichCells(object = subset(OB_s,subset = sample =="OB1"& center.x < xmax & center.x > xmin & center.y > ymin & center.y < ymax))],OB_s@active.ident[WhichCells(object = subset(OB_s,subset = sample =="OB1"& center.x < xmax & center.x > xmin & center.y > ymin & center.y < ymax))])
x <- OB_s$center.x[WhichCells(object = subset(OB_s,subset = sample =="OB1"& center.x < xmax & center.x > xmin & center.y > ymin & center.y < ymax))]
y <- OB_s$center.y[WhichCells(object = subset(OB_s,subset = sample =="OB1"& center.x < xmax & center.x > xmin & center.y > ymin & center.y < ymax))]
group <- OB_s@active.ident[WhichCells(object = subset(OB_s,subset = sample =="OB1"& center.x < xmax & center.x > xmin & center.y > ymin & center.y < ymax))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=0.2)+ theme_linedraw()+coord_fixed(ratio = 1)+theme(legend.key.size = unit(0.01, 'cm'))+theme(legend.position="bottom", legend.text = element_text(size=5))
#ggplot(dfall,aes(x,y,colour=group))+geom_point(size=(2000/(xmax-xmin)))+ theme_linedraw()+coord_fixed(ratio = 1)+theme(legend.key.size = unit(0.01, 'cm'))+theme(legend.position="bottom", legend.text = element_text(size=5))
ggsave('file2.pdf')
#features=c('Gad1', 'Slc17a6', 'Ttyh2', 'Mbp', 'Pdgfra', 'Aqp4', 'Selplg', 'Cd24a', 'Fn1', 'Myh11')
features=c("Tmem119","P2ry12","Aif1","Cx3cr1","Ptprc","Cd68","C3","Csf1r","Itgax","Clec7a","Spp1","Lyve1","Cd163","Mrc1","Ccr2","Camk2a","Snap25","Gad1","Gad2","Cd79a","Cldn5","Pecam1","Fn1","Tek","Cd4","Cd3e","Cd8a","Itgae","Ighmbp2","Des","Pdgfrb","Cspg4","Pdgfra","Aspa","Plp1","Aldh1l1","Aqp4","Gfap")
DotPlot(OB_s, features = features) + RotatedAxis()
#plot XY, color by gene expression
#gene <- "	Olfr1271"
#sampleA <- "OB1"
#dfall<- data.frame(OB.combined$center.x[WhichCells(object = subset(OB.combined,subset = sample ==sampleA))],OB.combined$center.y[WhichCells(object = subset(OB.combined,subset = sample ==sampleA))],FetchData(object =subset(OB.combined,subset = sample ==sampleA),vars =gene,slot="counts"))
#x <- OB.combined$center.x[WhichCells(object = subset(OB.combined,subset = sample ==sampleA))]
#y <- OB.combined$center.y[WhichCells(object = subset(OB.combined,subset = sample ==sampleA))]
#Expression <-dfall[,3]
#mid<-mean(Expression)
#ggplot(dfall,aes(x,y,color=Expression))+geom_point(size=.1)+ggtitle(paste(sampleA,gene,sep =" "))+ theme_linedraw()+scale_color_gradient(trans="log")
