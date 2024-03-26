#Setup the Seurat Object
#Install pakages (only once)
install.packages('Seurat')
install.packages("devtools")
install.packages("umap-learn")
install.packages("xlsx")
install.packages("ggpubr")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dittoSeq")
install.packages("matrixStats")
install.packages("ggdendro")
devtools::install_github("YuLab-SMU/ggtree")
install.packages("patchwork")
install.packages('PupillometryR')
#Activate the following libraries/packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(devtools)
library(tidyverse)
library(xlsx)
library(ggpubr)
library(dittoSeq)
library(matrixStats)
library(ggdendro)
library(cowplot)
#library(ggtree)
library(patchwork) 
library(PupillometryR)
# Load the pilot dataset
pilot.data<-readRDS("F:/Marta/Getting started with scRNAseq analysis/Pilot results/scd-22aeac/diagnostics/2021-10-11-JEN-AU.rds",refhook = NULL)

#remove chromosome lables
row.names(pilot.data)<-gsub("__.....","",row.names(pilot.data))
row.names(pilot.data)<-gsub("__....","",row.names(pilot.data))

# in seurat we will not make use of the spike-ins, so remove them from the expression matrix before creating the Seurat object. 
## Calculate ERCC total reads and the abundances (%) on the raw counts before creating a Seurat object
ercc <- grep("ERCC-",rownames(pilot.data))
total_ercc_reads<-Matrix::colSums(pilot.data[ercc,])
percent.ERCC <- (Matrix::colSums(pilot.data[ercc,])/Matrix::colSums(pilot.data))*100
mito <- grep("^mt",rownames(pilot.data))
total_mito_reads<-Matrix::colSums(pilot.data[mito,])
pilot.data.matrix<-as.matrix(pilot.data)
Var_expression_per_cell<-colVars(pilot.data.matrix)
Mean_expression_per_cell<-colMeans(pilot.data)

#when we create the seurat object we need to ask to remove the transcripts belonging to the ercc matrix and 
#add total_ercc_reads and percent.ERCC to object@meta.data in the total_ercc_reads and percent.ERCC column respectively
pilot <- CreateSeuratObject(counts = pilot.data[-ercc,], project = "Pilot analysis", min.cells = 0, min.features = 0)
pilot


pilot$orig.ident <- factor(pilot$orig.ident, levels=c("JEN-AU-s004","JEN-AU-s002"))

current.ids <- c("JEN-AU-s004","JEN-AU-s002")
new.ids <- c("Fresh_cells","Fixed_cells")
pilot@active.ident <- plyr::mapvalues(x = pilot@active.ident, from = current.ids, to = new.ids)
pilot$Condition <- Idents(pilot)

pilot$Condition <- factor(x=pilot$Condition, levels=c("Fresh_cells","Fixed_cells"))


#Pre-processing workflow

#Preparing all the extra columns in our surat Object
#QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pilot[["Var_expression_per_cell"]] <- Var_expression_per_cell
pilot[["Mean_expression_per_cell"]] <- Mean_expression_per_cell
pilot[["log_Var_expression_per_cell"]] <- log10(pilot$Var_expression_per_cell)
pilot[["log_Mean_expression_per_cell"]] <- log10(pilot$Mean_expression_per_cell)
pilot[["total_mito_reads"]] <- total_mito_reads
pilot[["total_ercc_reads"]] <- total_ercc_reads
pilot[["percent_ERCC"]] <- percent.ERCC
pilot[["percent_mt"]] <- PercentageFeatureSet(pilot, pattern = "^mt")
pilot[["percent_cyto"]] <- 100-(pilot$percent_mt)
pilot[["percent_ERCC"]] <- PercentageFeatureSet(pilot, pattern = "^ERCC")
pilot[["log10_nCount_RNA"]] <- log10(pilot$nCount_RNA)
pilot[["log10_nFeature_RNA"]] <- log10(pilot$nFeature_RNA)
pilot[["percent_ERCC_after_removal"]] <- PercentageFeatureSet(pilot, pattern = "^ERCC")
pilot[["ratio_mol_gen"]]<-pilot$nCount_RNA/pilot$nFeature_RNA
Control_of_ERCC_removal<- summary(pilot$percent_ERCC_after_removal)
Control_of_ERCC_removal
Control_percent_mt<-summary(pilot$percent_mt)
Control_percent_mt



pilot.list<-SplitObject(pilot, split.by = "Condition")

Fresh_cells<-pilot.list$Fresh_cells
Fixed_cells<-pilot.list$Fixed_cells

#Included and excluded labelling
Included_cells<-ifelse(Fresh_cells$nCount_RNA > 800 & Fresh_cells$nCount_RNA < 35000 & Fresh_cells$total_ercc_reads >500 &
                           Fresh_cells$nFeature_RNA >500  & Fresh_cells$ratio_mol_gen>1.2 & Fresh_cells$percent_mt<100,"INcluded","EXcluded")
Fresh_cells[["Included_cells"]]<-Included_cells

Included_cells<-ifelse(Fixed_cells$nCount_RNA > 800 & Fixed_cells$nCount_RNA < 35000 & Fixed_cells$total_ercc_reads >500 &
                         Fixed_cells$nFeature_RNA >500  & Fixed_cells$ratio_mol_gen>1.2 & Fixed_cells$percent_mt<100,"INcluded","EXcluded")
Fixed_cells[["Included_cells"]]<-Included_cells


#Scatter plots to visualise - Included_cells (ALL BARCODES)

#---Fresh
plot1 <- FeatureScatter(Fresh_cells, feature1 = "log10_nCount_RNA", feature2 = "percent_mt", group.by = "Included_cells", plot.cor=F) + xlim(0, 5) + ylim(0, 100)+ geom_vline(xintercept = 2.92942) + geom_vline(xintercept = 4.54407) + geom_hline(yintercept=100)+ theme(legend.position="none")
plot2 <- FeatureScatter(Fresh_cells, feature1 = "log10_nFeature_RNA", feature2 = "percent_mt", group.by = "Included_cells", plot.cor=F) + xlim(0, 4) + ylim(0, 100)+geom_vline(xintercept = 2.69897) + theme(legend.position="none")
plot3 <- FeatureScatter(Fresh_cells, feature1 = "log10_nFeature_RNA", feature2 = "log10_nCount_RNA", group.by = "Included_cells", plot.cor=F) + xlim(0, 4) + ylim(0, 5)+ geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot4 <- FeatureScatter(Fresh_cells, feature1 = "ratio_mol_gen", feature2 = "log10_nCount_RNA", group.by = "Included_cells", plot.cor=F) + ylim(0, 5) + geom_vline(xintercept = 1.2) + geom_hline(yintercept=2.92942)+ geom_hline(yintercept = 4.54407)  + theme(legend.position="none")
plot5 <- FeatureScatter(Fresh_cells, feature1 = "ratio_mol_gen", feature2 = "log10_nFeature_RNA", group.by = "Included_cells", plot.cor=F) + ylim(0, 4)+ geom_vline(xintercept = 1.2) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot6 <- FeatureScatter(Fresh_cells, feature1 = "total_ercc_reads", feature2 = "log10_nFeature_RNA", group.by = "Included_cells", plot.cor=F) + xlim(0, 2500) + ylim(0, 4) + geom_vline(xintercept = 500) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot7 <- FeatureScatter(Fresh_cells, feature1 = "total_ercc_reads", feature2 = "log10_nCount_RNA", group.by = "Included_cells", plot.cor=F) + xlim(0, 2500) + ylim(0, 5) + geom_vline(xintercept = 500) + geom_hline(yintercept=2.92942)+ geom_hline(yintercept = 4.54407)  + theme(legend.position="none")

annotate_figure(ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7), top = text_grob("Fresh"))

#---Fixed
plot1 <- FeatureScatter(Fixed_cells, feature1 = "log10_nCount_RNA", feature2 = "percent_mt", group.by = "Included_cells", plot.cor=F) + xlim(0, 5) + ylim(0, 100)+ geom_vline(xintercept = 2.92942) + geom_vline(xintercept = 4.54407) + geom_hline(yintercept=100)+ theme(legend.position="none")
plot2 <- FeatureScatter(Fixed_cells, feature1 = "log10_nFeature_RNA", feature2 = "percent_mt", group.by = "Included_cells", plot.cor=F) + xlim(0, 4) + ylim(0, 100)+geom_vline(xintercept = 2.69897) + theme(legend.position="none")
plot3 <- FeatureScatter(Fixed_cells, feature1 = "log10_nFeature_RNA", feature2 = "log10_nCount_RNA", group.by = "Included_cells", plot.cor=F) + xlim(0, 4) + ylim(0, 5)+ geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot4 <- FeatureScatter(Fixed_cells, feature1 = "ratio_mol_gen", feature2 = "log10_nCount_RNA", group.by = "Included_cells", plot.cor=F) + ylim(0, 5) + geom_vline(xintercept = 1.2) + geom_hline(yintercept=2.92942)+ geom_hline(yintercept = 4.54407)  + theme(legend.position="none")
plot5 <- FeatureScatter(Fixed_cells, feature1 = "ratio_mol_gen", feature2 = "log10_nFeature_RNA", group.by = "Included_cells", plot.cor=F) + ylim(0, 4)+ geom_vline(xintercept = 1.2) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot6 <- FeatureScatter(Fixed_cells, feature1 = "total_ercc_reads", feature2 = "log10_nFeature_RNA", group.by = "Included_cells", plot.cor=F) + xlim(0, 2500) + ylim(0, 4) + geom_vline(xintercept = 500) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot7 <- FeatureScatter(Fixed_cells, feature1 = "total_ercc_reads", feature2 = "log10_nCount_RNA", group.by = "Included_cells", plot.cor=F) + xlim(0, 2500) + ylim(0, 5) + geom_vline(xintercept = 500) + geom_hline(yintercept=2.92942)+ geom_hline(yintercept = 4.54407)  + theme(legend.position="none")

annotate_figure(ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7), top = text_grob("Fixed"))

#Scatter plots to visualise - Included_cells (ALL BARCODES)
Fresh_meta_df<-Fresh_cells@meta.data

ggplot(Fresh_meta_df, aes(x=log10_nFeature_RNA, y=log10_nCount_RNA, color=percent_mt)) + geom_point()+ xlim(0, 4) + ylim(0, 5)+
  scale_color_gradient(low="gray", high="blue", limits = c(0, 100))+ theme_classic() + geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897)+
  ggtitle("Fresh")

ggplot(Fresh_meta_df, aes(x=log10_nFeature_RNA, y=log10_nCount_RNA, color=ratio_mol_gen)) + geom_point()+ xlim(0, 4) + ylim(0, 5)+
  scale_color_gradient(low="gray", high="blue")+ theme_classic() + geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897)+
  ggtitle("Fresh")

ggplot(Fresh_meta_df, aes(x=log10_nFeature_RNA, y=log10_nCount_RNA, color=total_mito_reads)) + geom_point()+ xlim(0, 4) + ylim(0, 5)+
  scale_color_gradient(low="gray", high="blue", limits = c(0, 3789))+ theme_classic() + geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897)+
  ggtitle("Fresh")

Fixed_meta_df<-Fixed_cells@meta.data

ggplot(Fixed_meta_df, aes(x=log10_nFeature_RNA, y=log10_nCount_RNA, color=percent_mt)) + geom_point()+ xlim(0, 4) + ylim(0, 5)+
  scale_color_gradient(low="gray", high="blue", limits = c(0, 100))+ theme_classic() + geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897)+
  ggtitle("Fixed")

ggplot(Fixed_meta_df, aes(x=log10_nFeature_RNA, y=log10_nCount_RNA, color=ratio_mol_gen)) + geom_point()+ xlim(0, 4) + ylim(0, 5)+
  scale_color_gradient(low="gray", high="blue")+ theme_classic() + geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897)+
  ggtitle("Fixed")

ggplot(Fixed_meta_df, aes(x=log10_nFeature_RNA, y=log10_nCount_RNA, color=total_mito_reads)) + geom_point()+ xlim(0, 4) + ylim(0, 5)+
  scale_color_gradient(low="gray", high="blue", limits = c(0, 3789))+ theme_classic() + geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897)+
  ggtitle("Fixed")

#Labeling included and excluded only with nFeatures threshold
#Included and excluded labelling
Included_cells_F<-ifelse(Fresh_cells$nFeature_RNA >500,"INcluded","EXcluded")
Fresh_cells[["Included_cells_F"]]<-Included_cells_F

Included_cells_F<-ifelse(Fixed_cells$nFeature_RNA >500,"INcluded","EXcluded")
Fixed_cells[["Included_cells_F"]]<-Included_cells_F

table(Fresh_cells@meta.data[["Included_cells"]])
table(Fixed_cells@meta.data[["Included_cells"]])
table(Fresh_cells@meta.data[["Included_cells_F"]])
table(Fixed_cells@meta.data[["Included_cells_F"]])

#---Fresh
plot1 <- FeatureScatter(Fresh_cells, feature1 = "log10_nCount_RNA", feature2 = "percent_mt", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 5) + ylim(0, 100)+ geom_vline(xintercept = 2.92942) + geom_vline(xintercept = 4.54407) + geom_hline(yintercept=100)+ theme(legend.position="none")
plot2 <- FeatureScatter(Fresh_cells, feature1 = "log10_nFeature_RNA", feature2 = "percent_mt", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 4) + ylim(0, 100)+geom_vline(xintercept = 2.69897) + theme(legend.position="none")
plot3 <- FeatureScatter(Fresh_cells, feature1 = "log10_nFeature_RNA", feature2 = "log10_nCount_RNA", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 4) + ylim(0, 5)+ geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot4 <- FeatureScatter(Fresh_cells, feature1 = "ratio_mol_gen", feature2 = "log10_nCount_RNA", group.by = "Included_cells_F", plot.cor=F) + ylim(0, 5) + geom_vline(xintercept = 1.2) + geom_hline(yintercept=2.92942)+ geom_hline(yintercept = 4.54407)  + theme(legend.position="none")
plot5 <- FeatureScatter(Fresh_cells, feature1 = "ratio_mol_gen", feature2 = "log10_nFeature_RNA", group.by = "Included_cells_F", plot.cor=F) + ylim(0, 4)+ geom_vline(xintercept = 1.2) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot6 <- FeatureScatter(Fresh_cells, feature1 = "total_ercc_reads", feature2 = "log10_nFeature_RNA", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 2500) + ylim(0, 4) + geom_vline(xintercept = 500) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot7 <- FeatureScatter(Fresh_cells, feature1 = "total_ercc_reads", feature2 = "log10_nCount_RNA", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 2500) + ylim(0, 5) + geom_vline(xintercept = 500) + geom_hline(yintercept=2.92942)+ geom_hline(yintercept = 4.54407)  + theme(legend.position="none")

annotate_figure(ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7), top = text_grob("Fresh"))

#---Fixed
plot1 <- FeatureScatter(Fixed_cells, feature1 = "log10_nCount_RNA", feature2 = "percent_mt", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 5) + ylim(0, 100)+ geom_vline(xintercept = 2.92942) + geom_vline(xintercept = 4.54407) + geom_hline(yintercept=100)+ theme(legend.position="none")
plot2 <- FeatureScatter(Fixed_cells, feature1 = "log10_nFeature_RNA", feature2 = "percent_mt", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 4) + ylim(0, 100)+geom_vline(xintercept = 2.69897) + theme(legend.position="none")
plot3 <- FeatureScatter(Fixed_cells, feature1 = "log10_nFeature_RNA", feature2 = "log10_nCount_RNA", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 4) + ylim(0, 5)+ geom_vline(xintercept = 2.69897)+ geom_hline(yintercept= 4.54407) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot4 <- FeatureScatter(Fixed_cells, feature1 = "ratio_mol_gen", feature2 = "log10_nCount_RNA", group.by = "Included_cells_F", plot.cor=F) + ylim(0, 5) + geom_vline(xintercept = 1.2) + geom_hline(yintercept=2.92942)+ geom_hline(yintercept = 4.54407)  + theme(legend.position="none")
plot5 <- FeatureScatter(Fixed_cells, feature1 = "ratio_mol_gen", feature2 = "log10_nFeature_RNA", group.by = "Included_cells_F", plot.cor=F) + ylim(0, 4)+ geom_vline(xintercept = 1.2) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot6 <- FeatureScatter(Fixed_cells, feature1 = "total_ercc_reads", feature2 = "log10_nFeature_RNA", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 2500) + ylim(0, 4) + geom_vline(xintercept = 500) + geom_hline(yintercept=2.69897) + theme(legend.position="none")
plot7 <- FeatureScatter(Fixed_cells, feature1 = "total_ercc_reads", feature2 = "log10_nCount_RNA", group.by = "Included_cells_F", plot.cor=F) + xlim(0, 2500) + ylim(0, 5) + geom_vline(xintercept = 500) + geom_hline(yintercept=2.92942)+ geom_hline(yintercept = 4.54407)  + theme(legend.position="none")

annotate_figure(ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7), top = text_grob("Fixed"))
