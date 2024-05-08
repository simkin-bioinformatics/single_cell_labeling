####USER INPUT 4####
####Differential Expression####

library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)

STRESS.combined = readRDS(file=snakemake@input[[1]])

STRESS.combined$celltype.Sample2<-paste(Idents(STRESS.combined), STRESS.combined$Sample2, sep="_")
STRESS.combined$celltype<-Idents(STRESS.combined)
Idents(STRESS.combined)<-"celltype.Sample2"


Neuron_stress<-FindMarkers(STRESS.combined, ident.1 = snakemake@params[[1]], ident.2 = snakemake@params[[2]], verbose=FALSE)
write.csv(Neuron_stress, file = snakemake@output[[1]])


neuron_data<-read.csv(snakemake@output[[1]], stringsAsFactors = TRUE)  
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

write.csv(neuron_data, file = snakemake@output[[1]])

#neuron_data<-read.csv("15_T5_vs_T0.csv", stringsAsFactors = TRUE)

library(ggrepel)
p5<-ggplot(data= neuron_data, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel))+
  geom_point()+
  theme_minimal()+
  ggrepel::geom_text_repel()+
  scale_color_manual(values=c("blue", "black", "red"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="red")+
  geom_hline(yintercept = -log10(0.05), col="red")

pdf(file=snakemake@output[[2]], width = 16, height = 16)
p5 + labs(title= snakemake@params[3])+ theme(plot.title = element_text(hjust=0.5))
dev.off()
#print(p5)
write.csv(neuron_data, file = snakemake@output[[1]])


