#Setup the Seurat Object
#Install pakages (only once)
install.packages('Seurat')
install.packages("devtools")
install.packages("umap-learn")
install.packages("xlsx")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dittoSeq")
install.packages("matrixStats")
install.packages("data.table")
install.packages("magrittr")
install.packages('clustree')
install.packages('egg')
#Activate the following libraries/packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(devtools)
library(tidyverse)
library(xlsx)
library(dittoSeq)
library(matrixStats)
library(data.table)
library(magrittr)
library(clustree)
library(egg)
library(pheatmap)
#-------------------------------------Load the two SO----------------------------------------------
#Hannah_A_SO_w_metadata
Hannah_Fresh<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/Hannah_A_SO_w_metadata.rds")
#P28WT_dec_filt2_doublets_to_merge_SO
P28WT_Fixed_with_doublets<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/P28WT_dec_filt2_doublets_to_merge_SO.rds")
P28WT_Fixed_with_doublets
#------------------------------------Prepare the objects------------------------------

#P24WT_Fresh

#Extract P24 from Hannah_A_SO_w_metadata
P24WT_Fresh<-subset(Hannah_Fresh, subset = Age == "24"| Age =="24*")
#Adding condition in metadata
P24WT_Fresh[["Condition"]] <- rep("Fresh",each=1063)
P24WT_Fresh[["percent_ribo"]] <- PercentageFeatureSet(P24WT_Fresh, pattern = "^Rp")
Idents(object = P24WT_Fresh) <- "Condition"
#Save
saveRDS(P24WT_Fresh,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/Hannah_A_SO_w_metadata.rds")
#Load
P24WT_Fresh<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/Hannah_A_SO_w_metadata.rds")


#P28WT_Fixed

#Clean metadata from P28WT_dec_filt2_doublets_to_merge_SO
#Create Seurat Object with RNA assay from the previous P28
RNA_P28<-GetAssay(object = P28WT_Fixed_with_doublets, assay = "RNA")
P28WT_Fixed <- CreateSeuratObject(counts = RNA_P28, project = "methanol")
P28WT_Fixed
#Add the rest of the assays (non-decontaminated)
#P28WT_Fixed[["RNA_no_dec"]]<-GetAssay(object = P28WT_Fixed_with_doublets, assay = "RNA_no_dec")

#Export the metadata as a data frame
meta_data_P28<-P28WT_Fixed_with_doublets@meta.data
colnames(meta_data_P28)

#Create new metadata
P28WT_Fixed[["Age"]] <- meta_data_P28$orig.ident
P28WT_Fixed[["percent_mt"]] <- PercentageFeatureSet(P28WT_Fixed, pattern = "^mt")
P28WT_Fixed[["total_mito_reads"]] <- meta_data_P28$total_mito_reads
P28WT_Fixed[["percent_ribo"]] <- PercentageFeatureSet(P28WT_Fixed, pattern = "^Rp")
P28WT_Fixed[["log10_nCount_RNA"]] <- log10(P28WT_Fixed$nCount_RNA)
P28WT_Fixed[["log10_nFeature_RNA"]] <- log10(P28WT_Fixed$nFeature_RNA)
P28WT_Fixed[["ratio_mol_gen"]]<-P28WT_Fixed$nCount_RNA/P28WT_Fixed$nFeature_RNA
#P28WT_Fixed[["X1st_cont_fract"]]<-meta_data_P28$X1st_cont_fract
#P28WT_Fixed[["X2nd_cont_fract"]]<-meta_data_P28$X2nd_cont_fract
#P28WT_Fixed[["DoubletFinder"]]<-meta_data_P28$DoubletFinder
P28WT_Fixed[["scDblFinder"]]<-meta_data_P28$scDblFinder
P28WT_Fixed[["Condition"]] <- rep("Fixed",each=5501)

DefaultAssay(P28WT_Fixed) <- "RNA"
Idents(object = P28WT_Fixed) <- "Condition"

#Save seurat object with clean metadata
saveRDS(P28WT_Fixed,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/P28WT_Fixed_metadata_clean.rds")
#Load seurat object with clean metadata
#P28WT_Fixed<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/P28WT_Fixed_metadata_clean.rds")
table(P28WT_Fixed$scDblFinder)
#Removing doublets
P28WT_Fixed<-subset(P28WT_Fixed, subset = scDblFinder=="singlet" & P28WT_Fixed$nCount_RNA < 35000)
P28WT_Fixed[["scDblFinder"]]<-NULL
P28WT_Fixed
#Save seurat object without doublets
saveRDS(P28WT_Fixed,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/P28WT_Fixed_metadata_clean_wo_doublets.rds")
#Load seurat object without doublets
#P28WT_Fixed<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/P28WT_Fixed_metadata_clean_wo_doublets.rds")
print("done")

#--------------------------------------------------------Normalisation----------------------------------------------------------------------------
P24WT_Fresh<-SCTransform(P24WT_Fresh,new.assay.name="SCT", variable.features.n=2000)
P24WT_Fresh

P28WT_Fixed<-SCTransform(P28WT_Fixed,new.assay.name="SCT", variable.features.n=2000)
P28WT_Fixed

#Save seurat object normalised
saveRDS(P24WT_Fresh,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/P24WT_Fresh_norm.rds")
saveRDS(P28WT_Fixed,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/P28WT_Fixed_norm.rds")
#Load seurat object normalised
#P24WT_Fresh<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/P24WT_Fresh_norm.rds")
P28WT_Fixed<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/P28WT_Fixed_norm.rds")
print("done")
#---------------------------------------------------------Integration --------------------------------------------------------------
options(future.globals.maxSize = 16000 * 1024^2)
#Create a list of Seurat objects 
methanol.list<-list(P24WT_Fresh,P28WT_Fixed)
# select features for downstream integration
methanol.features <- SelectIntegrationFeatures(object.list = methanol.list, nfeatures = 2000)
#run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated
methanol.list<-PrepSCTIntegration(object.list = methanol.list, anchor.features = methanol.features, 
                                verbose = FALSE)
#identify anchors
methanol.anchors<-FindIntegrationAnchors(object.list = methanol.list, normalization.method = "SCT", 
                                      anchor.features = methanol.features, verbose = FALSE)
#integrate the datasets
methanol.integrated<-IntegrateData(anchorset = methanol.anchors, normalization.method = "SCT", 
                                verbose = FALSE)
#Save methanol integrated
saveRDS(methanol.integrated,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/methanol.integrated.rds")
print("done")

#once integrated we need to re-assign the order of the conditions
methanol.integrated$Condition <- factor(methanol.integrated$Condition, levels=c("Fresh","Fixed"))
Idents(methanol.integrated) <- "Condition"

#--------------------------------------------------Perform linear dimensional reduction-------------------------------------------------
methanol.integrated <- RunPCA(methanol.integrated, features = VariableFeatures(object = methanol.integrated))
# Examine and visualize PCA results a few different ways
DimPlot(methanol.integrated, reduction = "pca", group.by = "Condition")
#Determine the ‘dimensionality’ of the dataset
ElbowPlot(methanol.integrated)

#--------------------------------------------------UMAP--------------------------------------------------
methanol.integrated <- RunUMAP(methanol.integrated, dims = 1:15, n.neighbors = 15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(methanol.integrated, reduction = "umap", label=FALSE)


#Save methanol integrated
saveRDS(methanol.integrated,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/methanol.integrated_UMAP.rds")

#Load methanol integrated
methanol.integrated <- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/methanol.integrated_UMAP.rds")


#-----------------------------------------Clustering and cluster inspection

#         1. Represent our cells on UMAP representation + plot different QC parameters

#UMAP
methanol.integrated <- RunUMAP(methanol.integrated, dims = 1:15, n.neighbors = 5, reduction.name= "UMAP_n5", reduction = "pca")
methanol.integrated <- RunUMAP(methanol.integrated, dims = 1:15, n.neighbors = 10, reduction.name= "UMAP_n10", reduction = "pca")
methanol.integrated <- RunUMAP(methanol.integrated, dims = 1:15, n.neighbors = 50, reduction.name= "UMAP_n50", reduction = "pca")
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
plot1<-DimPlot(methanol.integrated, reduction = "UMAP_n5", label=T, pt.size=1) + ggtitle("UMAP_n5")
plot2<-DimPlot(methanol.integrated, reduction = "UMAP_n10", label=T, pt.size=1)+ ggtitle("UMAP_n10")
plot3<-DimPlot(methanol.integrated, reduction = "UMAP_n50", label=T, pt.size=1)+ ggtitle("UMAP_n50")

ggarrange(plot1, plot2, plot3, nrow=2)

DimPlot(methanol.integrated, reduction = "UMAP_n10", label=T, pt.size=1)
VlnPlot(methanol.integrated,group.by = "Condition", features= c("nCount_RNA","nFeature_RNA", "X1st_cont_fract",
                                                                "log10_nCount_RNA", "log10_nFeature_RNA", "X2nd_cont_fract", 
                                                                "total_mito_reads", "percent_mt", "percent_ribo",
                                                                "ratio_mol_gen"), ncol = 4)

FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("nCount_RNA","nFeature_RNA", "X1st_cont_fract",
                                                                    "log10_nCount_RNA", "log10_nFeature_RNA", "X2nd_cont_fract", 
                                                                    "total_mito_reads", "percent_mt", "percent_ribo",
                                                                    "ratio_mol_gen"), ncol = 4)
FeaturePlot(methanol.integrated, reduction="pca", features= c("nCount_RNA","nFeature_RNA", "X1st_cont_fract",
                                                              "log10_nCount_RNA", "log10_nFeature_RNA", "X2nd_cont_fract", 
                                                              "total_mito_reads", "percent_mt", "percent_ribo",
                                                              "ratio_mol_gen"), ncol = 4)


#         2. Use known cell identity markers to see where are the cell identities represented in the visualization.


#To be able to see the markers and not have a to wide range of expression we need to plot the log2(1+counts) expression
#(if we had 1 cell with an expression of 40 we were no able to see cells with expression of 15)

#Create a new assay to store log2(1+counts) information
#Save the linear normalized (LN) counts matrix
LN_counts_methanol.integrated<-as.data.frame(methanol.integrated@assays[["SCT"]]@counts)
#add +1 to each element of the matrix
LN_counts_methanol.integrated_1<-LN_counts_methanol.integrated+1
#check that it worked
LN_counts_methanol.integrated[1:10,1:4]
LN_counts_methanol.integrated_1[1:10,1:4]
#do the log2 from the counts+1 matrix
LN_counts_methanol.integrated_1_log2<-log2(LN_counts_methanol.integrated_1)
LN_counts_methanol.integrated_1_log2[1:10,1:4]
#Create a new assay to store log2(1+counts) information
log2_assay <- CreateAssayObject(counts = LN_counts_methanol.integrated_1_log2)

# add this assay to the previously created Seurat object
methanol.integrated[["LOG2"]] <- log2_assay

# Validate that the object now contains multiple assays
Assays(methanol.integrated)
#Check the active assay 
DefaultAssay(methanol.integrated)
# Switch the default to LOG2
DefaultAssay(methanol.integrated) <- "LOG2"
DefaultAssay(methanol.integrated)

#Save seurat object with log2(normalised counts+1) 
saveRDS(methanol.integrated,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/methanol.integrated_log2.rds")

#Load seurat object with log2(normalised counts+1)
methanol.integrated<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/methanol.integrated_log2.rds")

#Now we can plot!
#(PILOT markers)
#Granule mature 
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Calb1","Plk5","Rasgrp1", "Smad3", "Kcnip3", "Mfsd4","Tanc1"), slot ="counts") + plot_annotation(title ="Granule mature")
#Granule immature
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Dsp", "Rasgrf2", "Gda", "Rbm24"), slot ="counts")+ plot_annotation(title ="Granule immature")
#Neuroblast 2
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Prox1","Dcx", "Sox11","Neurod1","Calb2"), slot ="counts")+ plot_annotation(title ="Neuroblast 2")
#Neuroblast1
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Prox1","Dcx", "Sox11", "Neurod1", "Igfbpl1", "Calb2","Eomes"), slot ="counts")+ plot_annotation(title ="Neuroblast1")
#nIPC 
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Prox1","Dcx","Neurod1", "Igfbpl1", "Ccnd2","Eomes","Tfap2c"), slot ="counts")+ plot_annotation(title ="nIPC")
#RGL
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Sox9","Sox2","Aldoc", "Id4", "Hopx", "Lpar1"), slot ="counts")+ plot_annotation(title ="RGL")
#RG
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Fam210b", "Vim", "Acot1","Notch2"), slot ="counts")+ plot_annotation(title ="RG (Notch2 as -RG and +RGL)")
#Astrocytes
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("S100b","Sox9","Sox2","Id3","Aqp4"), slot ="counts")+ plot_annotation(title ="Astrocytes")
#Cajal-Retzius
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Reln", "Lhx1", "Lhx5", "NhIh2", "TrP283", "Diablo", "Ndnf", "Cacna2d2"), slot ="counts")+ plot_annotation(title ="Cajal-Retzius")
#OPC 
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Olig1","Olig2","Sox10","Pdgfra","Cspg4"), slot ="counts")+ plot_annotation(title ="OPC")
#Oligo
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Lpar1","Mbp","Plp1"), slot ="counts")+ plot_annotation(title ="Oligo")
#Endothelial cells
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Prom1","Esam","Hspb1"), slot ="counts")+ plot_annotation(title ="Endothelial cells")
#Microglia
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("Csf1r" , "Cx3cr1","P2ry12", "Fcrls", "Tmem119", "Olfml3", "Hexb"), slot ="counts")+ plot_annotation(title ="Microglia")


DimPlot(methanol.integrated, reduction = "UMAP_n10", label=FALSE, pt.size=1, group.by="orig.ident") + ggtitle("UMAP_n10")

#         3.Play around with the K of the function find neighbors until we reach our desired resolution


DefaultAssay(methanol.integrated) <- "integrated"
DefaultAssay(methanol.integrated)
#Cluster the cells (8 dims) and K=50 and resolution =0.5. This gives us a visual reference to interpret easier the cluster trees 
methanol.integrated <- FindNeighbors(methanol.integrated, dims = 1:15, k.param=50 , reduction = "pca")
methanol.integrated<- FindClusters(methanol.integrated, resolution = 0.5, algorithm=1)
#UMAP
#methanol.integrated <- RunUMAP(methanol.integrated, dims = 1:15, n.neighbors = 10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(methanol.integrated, reduction = "UMAP_n10", label=T, pt.size=1)+ plot_annotation(title ="dim_20 K_50 res_0.5")

#Playing with resolution
methanol.integrated <- FindNeighbors(methanol.integrated, dims = 1:15, k.param=50, reduction = "pca")
methanol.integrated<- FindClusters(methanol.integrated, resolution = seq(from = 0.1, to = 1, by = 0.1), algorithm=1)

clustree(methanol.integrated, prefix = "integrated_snn_res.")

#Playing with K
methanol.integrated_K10 <- FindNeighbors(methanol.integrated, dims = 1:15, k.param=10, reduction = "pca")
methanol.integrated_K25 <- FindNeighbors(methanol.integrated, dims = 1:15, k.param=25, reduction = "pca")
methanol.integrated_K50 <- FindNeighbors(methanol.integrated, dims = 1:15, k.param=50, reduction = "pca")
methanol.integrated_K75 <- FindNeighbors(methanol.integrated, dims = 1:15, k.param=75, reduction = "pca")
methanol.integrated_K100 <- FindNeighbors(methanol.integrated, dims = 1:15, k.param=100, reduction = "pca")

methanol.integrated_K10<- FindClusters(methanol.integrated_K10, resolution = 0.5, algorithm=1)
methanol.integrated_K25<- FindClusters(methanol.integrated_K25, resolution = 0.5, algorithm=1)
methanol.integrated_K50<- FindClusters(methanol.integrated_K50, resolution = 0.5, algorithm=1)
methanol.integrated_K75<- FindClusters(methanol.integrated_K75, resolution = 0.5, algorithm=1)
methanol.integrated_K100<- FindClusters(methanol.integrated_K100, resolution = 0.5, algorithm=1)



methanol.integrated[["K10"]] <- methanol.integrated_K10@meta.data[["seurat_clusters"]]
methanol.integrated[["K25"]] <- methanol.integrated_K25@meta.data[["seurat_clusters"]]
methanol.integrated[["K50"]] <- methanol.integrated_K50@meta.data[["seurat_clusters"]]
methanol.integrated[["K75"]] <- methanol.integrated_K75@meta.data[["seurat_clusters"]]
methanol.integrated[["K100"]] <- methanol.integrated_K100@meta.data[["seurat_clusters"]]

clustree(methanol.integrated, prefix = "K")

DimPlot(methanol.integrated, reduction = "UMAP_n10", label=T, pt.size=1, group.by="integrated_snn_res.0.2") + ggtitle("UMAP_n10_dim_20_K_50_res_0.2")

#Save seurat object with clustering
saveRDS(methanol.integrated,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/methanol.integrated_clust.rds")

#Load seurat object with clustering
methanol.integrated<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/methanol.integrated_clust.rds")


#Violin plots of QGs depending on cluster
#K=50 and resolution =0.2
VlnPlot(methanol.integrated, features= c("nCount_RNA","nFeature_RNA", "X1st_cont_fract",
                                         "log10_nCount_RNA", "log10_nFeature_RNA", "X2nd_cont_fract", 
                                         "total_mito_reads", "percent_mt", "percent_ribo",
                                         "ratio_mol_gen"),group.by="integrated_snn_res.0.2",ncol = 5)
#UMAP of QC parameters
FeaturePlot(methanol.integrated,reduction = "UMAP_n10", features= c("nCount_RNA","nFeature_RNA", "X1st_cont_fract",
                                                                    "log10_nCount_RNA", "log10_nFeature_RNA", "X2nd_cont_fract", 
                                                                    "total_mito_reads", "percent_mt", "percent_ribo",
                                                                    "ratio_mol_gen"), ncol = 4)


Idents(object = methanol.integrated) <- "integrated_snn_res.0.2"

#Lists of differential expressed genes (clusters)
methanol.markers <- FindAllMarkers(object=methanol.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox",assay="RNA")

#FDR correction
methanol.markers$BH_p_val_adj<-p.adjust(methanol.markers$p_val,method="BH")
#cell_expression_ratio
methanol.markers$cell_expression_ratio<-methanol.markers$pct.1/methanol.markers$pct.2


time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(methanol.markers, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/", paste(time.name), "methanol.markers dim_20_K_50_res_0.2.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)


#Cluster identification with dim_20_K_50_res_0.3


DimPlot(methanol.integrated, reduction = "UMAP_n10", label=T, pt.size=1, group.by="integrated_snn_res.0.3") + ggtitle("UMAP_n10_dim_20_K_50_res_0.3")

VlnPlot(methanol.integrated, features= c("nCount_RNA","nFeature_RNA", "X1st_cont_fract",
                                         "log10_nCount_RNA", "log10_nFeature_RNA", "X2nd_cont_fract", 
                                         "total_mito_reads", "percent_mt", "percent_ribo",
                                         "ratio_mol_gen"),group.by="integrated_snn_res.0.3",ncol = 5)

Idents(object = methanol.integrated) <- "integrated_snn_res.0.3"

#Lists of differential expressed genes (clusters)
methanol.markers_0.3 <- FindAllMarkers(object=methanol.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox",assay="RNA")

#FDR correction
methanol.markers_0.3$BH_p_val_adj<-p.adjust(methanol.markers_0.3$p_val,method="BH")
#cell_expression_ratio
methanol.markers_0.3$cell_expression_ratio<-methanol.markers_0.3$pct.1/methanol.markers_0.3$pct.2


time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(methanol.markers_0.3, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/", paste(time.name), "methanol.markers dim_20_K_50_res_0.3.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)




DimPlot(methanol.integrated, reduction = "UMAP_n10", label=T, pt.size=1, group.by="Cluster") + ggtitle("Cluster_hannah")


#ONLY GRANULE NURONS

Granule_neurons.integrated<-subset(methanol.integrated, subset = Cluster==c("Granule-immature","Granule-mature", NA))


DimPlot(Granule_neurons.integrated, reduction = "UMAP_n10", label=T, pt.size=1, group.by="Cluster") + ggtitle("Cluster_hannah")

#After cluster annotation with Linnarsson's database

Idents(object = methanol.integrated) <- "integrated_snn_res.0.3"
#Assigning cell type identity to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
new.cluster.ids <- c("Granule neurons","Granule neurons","Granule neurons","Granule neurons","Astrocytes","nIPCs + Neuroblasts",
                     "Microglia", "Astrocytes +  RGL", "Mossy cells", "Endothelial", "RGL", "Oligo", "GABA", "OPC", "Cajal Retzius")
methanol.integrated@active.ident <- plyr::mapvalues(x = methanol.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
methanol.integrated$CellType <- Idents(methanol.integrated)


#Save seurat object with cell identities
saveRDS(methanol.integrated,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/methanol.integrated_cell_identities.rds")

#Load seurat object with cell identities
methanol.integrated<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/Seurat objects/methanol.integrated_cell_identities.rds")
table(methanol.integrated$Condition)

#stacked barplots to see Condition contribution to cluster formation
Idents(methanol.integrated) <- "Condition"
dittoBarPlot(methanol.integrated,var="Condition", group.by = "CellType",retain.factor.levels=TRUE, color.panel=c("#F8766D","#00BFC4"))


#pie graphs with % of the different cell identities
## extract meta data
md <- methanol.integrated@meta.data %>% as.data.table
# the resulting md object has one "row" per cell

## count the number of cells per unique combinations of "Condition" and "CellType" ad plot the pie chart
Proportion_CellType_Condition<-md[, .N, by = c("Condition", "CellType")]
Proportion_Fresh<-subset(Proportion_CellType_Condition, Condition=="Fresh")
bp_Fresh<- ggplot(Proportion_Fresh, aes(x="", y=N, fill=CellType))+
  geom_bar(width = 1, stat = "identity")
bp_Fresh
pie_Fresh <- bp_Fresh + coord_polar("y", start=0) + ggtitle("Proportion_Fresh")+ geom_text(aes(label=CellType),position = position_stack(vjust = 0.5), size=3.5)
pie_Fresh



Proportion_Fixed<-subset(Proportion_CellType_Condition, Condition=="Fixed")
bp_Fixed<- ggplot(Proportion_Fixed, aes(x="", y=N, fill=CellType))+
  geom_bar(width = 1, stat = "identity")
bp_Fixed
pie_Fixed <- bp_Fixed + coord_polar("y", start=0) + ggtitle("Proportion_Fixed") + geom_text(aes(label=CellType),position = position_stack(vjust = 0.5), size=3.5)
pie_Fixed


#Proportion of each cell type (%)
Proportion_CellType_Condition
#Calculate the percentage of each condition inside each cell type
Astrocytes_RGL_Prop<-subset(Proportion_CellType_Condition, CellType=="Astrocytes +  RGL")
Astrocytes_RGL_Prop$Percentage<-(Astrocytes_RGL_Prop$N/sum(Astrocytes_RGL_Prop$N))*100

RGL_Prop<-subset(Proportion_CellType_Condition, CellType=="RGL")
RGL_Prop$Percentage<-(RGL_Prop$N/sum(RGL_Prop$N))*100

Mossy_cells_Prop<-subset(Proportion_CellType_Condition, CellType=="Mossy cells")
Mossy_cells_Prop$Percentage<-(Mossy_cells_Prop$N/sum(Mossy_cells_Prop$N))*100

nIPCs_Neuroblasts_Prop<-subset(Proportion_CellType_Condition, CellType=="nIPCs + Neuroblasts")
nIPCs_Neuroblasts_Prop$Percentage<-(nIPCs_Neuroblasts_Prop$N/sum(nIPCs_Neuroblasts_Prop$N))*100

Granule_neurons_Prop<-subset(Proportion_CellType_Condition, CellType=="Granule neurons")
Granule_neurons_Prop$Percentage<-(Granule_neurons_Prop$N/sum(Granule_neurons_Prop$N))*100

Astrocytes_Prop<-subset(Proportion_CellType_Condition, CellType=="Astrocytes")
Astrocytes_Prop$Percentage<-(Astrocytes_Prop$N/sum(Astrocytes_Prop$N))*100

Microglia_Prop<-subset(Proportion_CellType_Condition, CellType=="Microglia")
Microglia_Prop$Percentage<-(Microglia_Prop$N/sum(Microglia_Prop$N))*100

Oligodendrocytes_Prop<-subset(Proportion_CellType_Condition, CellType=="Oligo")
Oligodendrocytes_Prop$Percentage<-(Oligodendrocytes_Prop$N/sum(Oligodendrocytes_Prop$N))*100

OPCs_Prop<-subset(Proportion_CellType_Condition, CellType=="OPCs")
OPCs_Prop$Percentage<-(OPCs_Prop$N/sum(OPCs_Prop$N))*100

Endothelial_Prop<-subset(Proportion_CellType_Condition, CellType=="Endothelial")
Endothelial_Prop$Percentage<-(Endothelial_Prop$N/sum(Endothelial_Prop$N))*100

Cajal_Retzius_Prop<-subset(Proportion_CellType_Condition, CellType=="Cajal Retzius")
Cajal_Retzius_Prop$Percentage<-(Cajal_Retzius_Prop$N/sum(Cajal_Retzius_Prop$N))*100

GABA_Prop<-subset(Proportion_CellType_Condition, CellType=="GABA")
GABA_Prop$Percentage<-(GABA_Prop$N/sum(GABA_Prop$N))*100

Proportion_CellType_Condition_Percentage<-rbind(Granule_neurons_Prop,Astrocytes_Prop,nIPCs_Neuroblasts_Prop,Microglia_Prop,Astrocytes_RGL_Prop,
                                                Mossy_cells_Prop,Endothelial_Prop, RGL_Prop, Oligodendrocytes_Prop, GABA_Prop, OPCs_Prop,Cajal_Retzius_Prop)

Proportion_CellType_Condition_Percentage

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Proportion_CellType_Condition_Percentage, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach//", paste(time.name), "Proportion_CellType_Condition with Percentage.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)



#Calculate the % of each celltype inside each condition

Fresh_Prop<-subset(Proportion_CellType_Condition, Condition=="Fresh")
Fresh_Prop$Percentage<-(Fresh_Prop$N/sum(Fresh_Prop$N))*100

Fixed_Prop<-subset(Proportion_CellType_Condition, Condition=="Fixed")
Fixed_Prop$Percentage<-(Fixed_Prop$N/sum(Fixed_Prop$N))*100

Proportion_inside_Condition_Percentage<-rbind(Fresh_Prop,Fixed_Prop)

Proportion_inside_Condition_Percentage

write.xlsx(Proportion_inside_Condition_Percentage, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach//Proportion_CellType_Condition with Percentage.xlsx"))
           , sheetName = "Sheet2", col.names = TRUE, row.names = TRUE, append = T)




#Calculating the pairwise correlation for all the cells we profiled but using the Pearson residuals

#Now I will do it with methanol.integrated@assays[["integrated"]]@scale.data

methanol_scaled_matrix<-as.matrix(methanol.integrated@assays[["integrated"]]@scale.data)

methanol_scaled_cor<-cor(methanol_scaled_matrix)
Condition_and_cluster_df<-data.frame(methanol.integrated$Condition)
Clusters<-as.factor(methanol.integrated$CellType)
Condition_and_cluster_df$Clusters<-Clusters
colnames(Condition_and_cluster_df)<-c("Condition", "Clusters")
my_colour<-list(Condition=c(Fresh = "#F8766D", Fixed = "#00BFC4"),
                Clusters=c(Astrocytes_RGL= "#F8766D", nIPCs_Neuroblasts= "#CD9600", Immature_neurons= "#7CAE00", 
                           Mature_neurons= "#00BE67",Microglia= "#00BFC4",Oligodendrocytes="#00A9FF", OPCs="#C77CFF", Endothelial="#FF61CC"))
pheatmap(methanol_scaled_cor,show_rownames=F,show_colnames=F,annotation_row =Condition_and_cluster_df,annotation_col =Condition_and_cluster_df,
         annotation_names_row = FALSE, annotation_names_col = FALSE,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")


annotation_colors=my_colour


#Subset cell identities
Astrocytes_RGL<-subset(methanol.integrated, subset = CellType=="Astrocytes +  RGL")
nIPCs_Neuroblasts<-subset(methanol.integrated, subset = CellType=="nIPCs + Neuroblasts")
Immature_neurons<-subset(methanol.integrated, subset = CellType=="Immature neurons")
Mature_neurons<-subset(methanol.integrated, subset = CellType=="Mature neurons")
Microglia<-subset(methanol.integrated, subset = CellType=="Microglia")
Oligodendrocytes<-subset(methanol.integrated, subset = CellType=="Oligodendrocytes")
OPCs<-subset(methanol.integrated, subset = CellType=="OPCs")
Endothelial_cells<-subset(methanol.integrated, subset = CellType=="Endothelial")

#number of cells per cell identity
table(Idents(pilot.integrated), pilot.integrated$Condition)

Astrocytes_RGL_matrix<-as.matrix(Astrocytes_RGL@assays[["RNA"]]@counts)

Astrocytes_RGL_cor<-cor(Astrocytes_RGL_matrix)
Condition_df_Astrocytes_RGL<-data.frame(Astrocytes_RGL$Condition)
colnames(Condition_df_Astrocytes_RGL)<-"Condition"
my_colour<-list(Condition=c(Fresh = "#F8766D", Fixed = "#00BFC4"))
pheatmap(Astrocytes_RGL_cor,show_rownames=F,show_colnames=F,annotation_row =Condition_df_Astrocytes_RGL,annotation_col =Condition_df_Astrocytes_RGL,
         annotation_colors=my_colour,annotation_names_row = FALSE, annotation_names_col = FALSE,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

nIPCs_Neuroblasts_matrix<-as.matrix(nIPCs_Neuroblasts@assays[["RNA"]]@counts)

nIPCs_Neuroblasts_cor<-cor(nIPCs_Neuroblasts_matrix)
Condition_df_nIPCs_Neuroblasts<-data.frame(nIPCs_Neuroblasts$Condition)
colnames(Condition_df_nIPCs_Neuroblasts)<-"Condition"
my_colour<-list(Condition=c(Fresh = "#F8766D", Fixed = "#00BFC4"))
pheatmap(nIPCs_Neuroblasts_cor,show_rownames=F,show_colnames=F,annotation_row =Condition_df_nIPCs_Neuroblasts,annotation_col =Condition_df_nIPCs_Neuroblasts,
         annotation_colors=my_colour,annotation_names_row = FALSE, annotation_names_col = FALSE,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")



#----------------------------------------------Correaltions tests------------------------------------------------------------

#log2(+1)=log2_1

methanol.bulk_RNA_counts_log2_1<-log2(1+(AverageExpression(methanol.integrated, group.by="Condition", slot= "counts"))$RNA)
methanol.bulk_RNA_counts_log2_1<-data.frame(methanol.bulk_RNA_counts_log2_1)

#Normal quantile plot

qqnorm(methanol.bulk_RNA_counts_log2_1$Fresh, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "methanol.bulk_RNA_counts_log2_1$Fresh_cells quantiles",na.rm=TRUE)
qqline(methanol.bulk_RNA_counts_log2_1$Fresh, col = "steelblue", lwd = 2)
qqnorm(methanol.bulk_RNA_counts_log2_1$Fixed, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "methanol.bulk_RNA_counts_log2_1$Fixed_cells quantiles")
qqline(methanol.bulk_RNA_counts_log2_1$Fixed, col = "steelblue", lwd = 2)

#Plot scatter plot

p1 <- ggplot(methanol.bulk_RNA_counts_log2_1, aes(Fixed,Fresh)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() +geom_smooth(method=lm, se =TRUE) + 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (log2(counts+1)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (log2(counts+1)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p1

p2 <- ggplot(methanol.bulk_RNA_counts_log2_1, aes(Fixed,Fresh)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() + geom_abline(intercept = 0, slope = 1, color = "blue")+ xlim(0, 9)+ylim(0, 9)+
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (log2(counts+1)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (log2(counts+1)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p2



#Run Spearman’s correlation test


cor.test(methanol.bulk_RNA_counts_log2_1$Fresh, methanol.bulk_RNA_counts_log2_1$Fixed, method=c("spearman"))
cor.test(methanol.bulk_RNA_counts_log2_1$Fresh, methanol.bulk_RNA_counts_log2_1$Fixed, method=c("pearson"))




#-------------------------------- PCA  plotting --------------------------------
Idents(methanol.integrated) <- "Condition"
P1<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 2), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                    axis.text=element_text(size=6),
                                                                                    axis.title.x = element_text(size =10),
                                                                                    axis.title.y = element_text(size =10))
P1
ggsave(
  filename="P1.pdf",
  plot = P1,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))

P2<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 3), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                    axis.text=element_text(size=6),
                                                                                    axis.title.x = element_text(size =10),
                                                                                    axis.title.y = element_text(size =10))
P2
ggsave(
  filename="P2.pdf",
  plot = P2,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P3<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 4), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                    axis.text=element_text(size=6),
                                                                                    axis.title.x = element_text(size =10),
                                                                                    axis.title.y = element_text(size =10))
P3
ggsave(
  filename="P3.pdf",
  plot = P3,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P4<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 5), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                    axis.text=element_text(size=6),
                                                                                    axis.title.x = element_text(size =10),
                                                                                    axis.title.y = element_text(size =10))
P4
ggsave(
  filename="P4.pdf",
  plot = P4,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P5<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 6), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                    axis.text=element_text(size=6),
                                                                                    axis.title.x = element_text(size =10),
                                                                                    axis.title.y = element_text(size =10))
P5
ggsave(
  filename="P5.pdf",
  plot = P5,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P6<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 7), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                    axis.text=element_text(size=6),
                                                                                    axis.title.x = element_text(size =10),
                                                                                    axis.title.y = element_text(size =10))
P6
ggsave(
  filename="P6.pdf",
  plot = P6,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P7<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 8), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                    axis.text=element_text(size=6),
                                                                                    axis.title.x = element_text(size =10),
                                                                                    axis.title.y = element_text(size =10))
P7
ggsave(
  filename="P7.pdf",
  plot = P7,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P8<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 9), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                    axis.text=element_text(size=6),
                                                                                    axis.title.x = element_text(size =10),
                                                                                    axis.title.y = element_text(size =10))
P8
ggsave(
  filename="P8.pdf",
  plot = P8,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P9<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 10), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                     axis.text=element_text(size=6),
                                                                                     axis.title.x = element_text(size =10),
                                                                                     axis.title.y = element_text(size =10))
P9
ggsave(
  filename="P9.pdf",
  plot = P9,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P10<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 11), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                      axis.text=element_text(size=6),
                                                                                      axis.title.x = element_text(size =10),
                                                                                      axis.title.y = element_text(size =10))
P10
ggsave(
  filename="P10.pdf",
  plot = P10,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P11<-DimPlot(methanol.integrated,order="Fresh",dims = c(3, 2), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                     axis.text=element_text(size=6),
                                                                                     axis.title.x = element_text(size =7),
                                                                                     axis.title.y = element_text(size =7))
P11
ggsave(
  filename="P11.pdf",
  plot = P11,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 1,
  width = 11.5,
  height = 10,
  units = c("cm"))
P12<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 12), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                                      axis.text=element_text(size=6),
                                                                                                      axis.title.x = element_text(size =7),
                                                                                                      axis.title.y = element_text(size =7))
P12

P13<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 13), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                                       axis.text=element_text(size=6),
                                                                                                       axis.title.x = element_text(size =7),
                                                                                                       axis.title.y = element_text(size =7))
P13
P14<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 14), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                                       axis.text=element_text(size=6),
                                                                                                       axis.title.x = element_text(size =7),
                                                                                                       axis.title.y = element_text(size =7))
P14
P15<-DimPlot(methanol.integrated,order="Fresh",dims = c(1, 15), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                                                       axis.text=element_text(size=6),
                                                                                                       axis.title.x = element_text(size =7),
                                                                                                       axis.title.y = element_text(size =7))
P15


ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10, ncol = 2, nrow = 5)




#UMAP graphs

UMAP_condition<-DimPlot(methanol.integrated, reduction = "UMAP_n10", label=F, pt.size=1, group.by="Condition")
ggsave(
  filename="UMAP_condition.pdf",
  plot = UMAP_condition,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 3,
  width = 11.5,
  height = 10,
  units = c("cm"))

UMAP_cellT<-DimPlot(methanol.integrated, reduction = "UMAP_n10", label=F, pt.size=1, group.by="CellType")
ggsave(
  filename="UMAP_cellT.pdf",
  plot = UMAP_cellT,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 3,
  width = 11.5,
  height = 10,
  units = c("cm"))

#pie graphs with % of the different cell identities
mdata <- as.data.frame(methanol.integrated@meta.data)
Cell_type<-mdata$CellType
mdata2<-as.data.frame(Cell_type)
mdata2$Condition<-mdata$Condition

Fresh_df<-subset(mdata2, Condition=="Fresh")
Fixed_df<-subset(mdata2, Condition=="Fixed")
#Fresh
Fresh_per_cent<-as.data.frame(table(Fresh_df$Cell_type))
Fresh_per_cent$per_cent<-Fresh_per_cent$Freq/sum(Fresh_per_cent$Freq)*100
Fresh_per_cent
Fr<-ggplot(Fresh_per_cent, aes(x="", y=per_cent, fill=Var1))+
  geom_bar(stat = "identity")+ coord_polar("y", start=0)+ theme_classic()+ggtitle("Fresh")

#Fixed
Fixed_per_cent<-as.data.frame(table(Fixed_df$Cell_type))
Fixed_per_cent$per_cent<-Fixed_per_cent$Freq/sum(Fixed_per_cent$Freq)*100
Fixed_per_cent
Fi<-ggplot(Fixed_per_cent, aes(x="", y=per_cent, fill=Var1))+
  geom_bar(stat = "identity")+ coord_polar("y", start=0)+ theme_classic()+ggtitle("Fixed")

Pie_graphs<-ggarrange(Fr,Fi)

ggsave(
  filename="Pie_graphs.pdf",
  plot = Pie_graphs,
  path = "F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/13 Hochgerner vs Urbach/R plots/",
  scale = 2,
  width = 11.5,
  height = 10,
  units = c("cm"))


#









































































#UMAP
methanol.integrated <- RunUMAP(methanol.integrated, dims = 1:15, n.neighbors = 15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(methanol.integrated, reduction = "umap", label=FALSE, pt.size=2.3)

#Differential expressed genes heatmap (clusters)
methanol.markers %>%
  group_by(cluster) %>%
  slice_max(n =10, order_by = avg_log2FC) -> top10_log2FC
top10_log2FC_HM<-DoHeatmap(methanol.integrated, features = top10_log2FC$gene, label=FALSE)
top10_log2FC_HM + labs(title="top10_log2FC
                       ")

#Differential expressed genes heatmap (clusters)
methanol.markers %>%
  group_by(cluster) %>%
  slice_min(n =10, order_by = p_val) -> top10_pval
top10_pval_HM<-DoHeatmap(methanol.integrated, features = top10_pval$gene, label=FALSE, assay = )
top10_pval_HM + labs(title="top10_pval
                       ")






#To see if methanol fixation maintain cell composition  Cell type composition graphs 
#(pie graphs with % of the different cell identities or graph of cluster contribution that we already have) (Denisenko 2020)













#UMAP plot where cells are colored by replicate?  
# First, store the current identities in a new column of meta.data called CellType
methanol.integrated$CellType <- Idents(methanol.integrated)
# Next, switch the identity class of all cells to reflect replicate ID
Idents(methanol.integrated) <- "Condition"
DimPlot(methanol.integrated, reduction = "umap", pt.size=2.3)
Idents(methanol.integrated) <- "CellType"
DimPlot(methanol.integrated, reduction = "umap", pt.size=2.3, split.by = "Condition")+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))




#stacked barplots to see Condition contribution to cluster formation
Idents(methanol.integrated) <- "Condition"
dittoBarPlot(methanol.integrated,var="Condition", group.by = "CellType",retain.factor.levels=TRUE, color.panel=c("#F8766D","#00BFC4"))

#pie graphs with % of the different cell identities
## extract meta data
md <- methanol.integrated@meta.data %>% as.data.table
# the resulting md object has one "row" per cell

## count the number of cells per unique combinations of "Condition" and "CellType" ad plot the pie chart
Proportion_CellType_Condition<-md[, .N, by = c("Condition", "CellType")]
Proportion_Fresh<-subset(Proportion_CellType_Condition, Condition=="Fresh")
bp_Fresh<- ggplot(Proportion_Fresh, aes(x="", y=N, fill=CellType))+
  geom_bar(width = 1, stat = "identity")
bp_Fresh
pie_Fresh <- bp_Fresh + coord_polar("y", start=0) + ggtitle("Proportion_Fresh")+ geom_text(aes(label=CellType),position = position_stack(vjust = 0.5), size=3.5)
pie_Fresh

Proportion_Fixed<-subset(Proportion_CellType_Condition, Condition=="Fixed")
bp_Fixed<- ggplot(Proportion_Fixed, aes(x="", y=N, fill=CellType))+
  geom_bar(width = 1, stat = "identity")
bp_Fixed
pie_Fixed <- bp_Fixed + coord_polar("y", start=0) + ggtitle("Proportion_Fixed") + geom_text(aes(label=CellType),position = position_stack(vjust = 0.5), size=3.5)
pie_Fixed

#Proportion of each cell type (%)
Proportion_CellType_Condition
#Calculate the percentage of each condition inside each cell type
Astrocytes_RGL_Prop<-subset(Proportion_CellType_Condition, CellType=="Astrocytes +  RGL")
Astrocytes_RGL_Prop$Percentage<-(Astrocytes_RGL_Prop$N/sum(Astrocytes_RGL_Prop$N))*100

nIPCs_Neuroblasts_Prop<-subset(Proportion_CellType_Condition, CellType=="nIPCs + Neuroblasts")
nIPCs_Neuroblasts_Prop$Percentage<-(nIPCs_Neuroblasts_Prop$N/sum(nIPCs_Neuroblasts_Prop$N))*100

Immature_neurons_Prop<-subset(Proportion_CellType_Condition, CellType=="Immature neurons")
Immature_neurons_Prop$Percentage<-(Immature_neurons_Prop$N/sum(Immature_neurons_Prop$N))*100

Mature_neurons_Prop<-subset(Proportion_CellType_Condition, CellType=="Mature neurons")
Mature_neurons_Prop$Percentage<-(Mature_neurons_Prop$N/sum(Mature_neurons_Prop$N))*100

Microglia_Prop<-subset(Proportion_CellType_Condition, CellType=="Microglia")
Microglia_Prop$Percentage<-(Microglia_Prop$N/sum(Microglia_Prop$N))*100

Oligodendrocytes_Prop<-subset(Proportion_CellType_Condition, CellType=="Oligodendrocytes")
Oligodendrocytes_Prop$Percentage<-(Oligodendrocytes_Prop$N/sum(Oligodendrocytes_Prop$N))*100

OPCs_Prop<-subset(Proportion_CellType_Condition, CellType=="OPCs")
OPCs_Prop$Percentage<-(OPCs_Prop$N/sum(OPCs_Prop$N))*100

Endothelial_Prop<-subset(Proportion_CellType_Condition, CellType=="Endothelial")
Endothelial_Prop$Percentage<-(Endothelial_Prop$N/sum(Endothelial_Prop$N))*100

Proportion_CellType_Condition_Percentage<-rbind(Astrocytes_RGL_Prop,nIPCs_Neuroblasts_Prop,Immature_neurons_Prop,
                                                Mature_neurons_Prop,Microglia_Prop,Oligodendrocytes_Prop,OPCs_Prop,Endothelial_Prop)

Proportion_CellType_Condition_Percentage

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Proportion_CellType_Condition_Percentage, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/methanol results/methanol analysis + results/9 See effect of fixation on cell identities/", paste(time.name), "Proportion_CellType_Condition with Percentage.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)


#Plots of with the typical cell type markers ONLY WITH FIXED CELLS!!!


#UMAP
methanol.integrated <- RunUMAP(methanol.integrated, dims = 1:15, n.neighbors = 15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(methanol.integrated, reduction = "umap", label=FALSE, pt.size=2.3)

methanol.list.int<-SplitObject(methanol.integrated, split.by = "Condition")

P28WT_Fixed.i<-methanol.list.int$P28WT_Fixed



DefaultAssay(P28WT_Fixed.i)<-"SCT"
#Granule mature 
FeaturePlot(P28WT_Fixed.i, features= c("Prox1","Plk5", "Ntng1", "Rbfox3", "Calb1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Prox1","Plk5", "Ntng1", "Rbfox3", "Calb1"), slot ="counts")
#Granule mature (round 2)
VlnPlot(P28WT_Fixed.i,ncol=6, features= c("Rgs4","Trpc6","Bc030499","Ntng1","Rasgrp1","Nos1","Snhg11","Smad3","Adarb2","Ipcef1","Kcnip3","Tekt5","Gm12216","Kart2","Hlf","Plekha2","Icam4","Mfsd4","Tanc1","Cstad"), slot ="counts")

#Granule immature
FeaturePlot(P28WT_Fixed.i, features= c("Ncam1","Prox1","Dcx", "Fxyd7","Rbfox3", "Calb2" ), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Ncam1","Prox1","Dcx", "Fxyd7","Rbfox3", "Calb2" ), slot ="counts")
#Granule immature (round 2)
VlnPlot(P28WT_Fixed.i,ncol=6, features= c("Camk4", "Stc1", "Sdk2", "Rasd1", "Fxyd7", "Il16", "Dsp", "Rspo3", "Prickle2", "Smim3", "Mcm6", "Rasgrf2", "Tenm1", "Scl4a4", "Gda", "Ptprj", "Rbm24","Omg", "Dnajc27"), slot ="counts")

#Neuroblast 2
FeaturePlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Gal", "Sox11","Neurod1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Gal", "Sox11","Neurod1"), slot ="counts")
#Neuroblast1
FeaturePlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Eomes", "Igfbpl1", "Calb2", "Neurog2", "Sox11", "Neurod4", "Tac2", "Neurod1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Eomes", "Igfbpl1", "Calb2", "Neurog2", "Sox11", "Neurod4", "Tac2", "Neurod1"), slot ="counts")
#nIPC 
FeaturePlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Ascl", "Cdk1", "Tfap2c", "Eomes", "Igfbpl1","Neurog2", "Top2a", "Neurod4", "Aurkb", "Neurod1", "Ccnd2"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Ascl", "Cdk1", "Tfap2c", "Eomes", "Igfbpl1","Neurog2", "Top2a", "Neurod4", "Aurkb", "Neurod1", "Ccnd2"), slot ="counts")
#RGL
FeaturePlot(P28WT_Fixed.i, features= c("Tfap2c", "Sox9", "Lpar1", "Ascl1","Sox2", "Nes", "Vnn1","Hes5", "Rhcg", "Scl1a3", "Wnt8b", "Aldoc", "Apoe", "Id4", "Hopx", "GFAP"), slot ="counts")
VlnPlot(P28WT_Fixed.i, features= c("Tfap2c", "Sox9", "Lpar1", "Ascl1","Sox2", "Nes", "Vnn1","Hes5", "Rhcg", "Scl1a3", "Wnt8b", "Aldoc", "Apoe", "Id4", "Hopx", "GFAP"), slot ="counts")
#Astrocytes
FeaturePlot(P28WT_Fixed.i, features= c("S100b", "Fzd2", "Sox9","Sox2","Hes5","Scl1a3"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("S100b", "Fzd2", "Sox9","Sox2","Hes5","Scl1a3"), slot ="counts")
#OPC 
FeaturePlot(P28WT_Fixed.i, features= c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"), slot ="counts")
#Oligo 
FeaturePlot(P28WT_Fixed.i, features= c("Lpar1","Mbp", "Plp1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Lpar1","Mbp", "Plp1"), slot ="counts")
#Endothelial cells
FeaturePlot(P28WT_Fixed.i, features= c("Nes", "Prom1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, ncol=6, features= c("Nes", "Prom1", "Efna1", "Tm4sf1", "Pcp4l1","Cd93","Igfbp7","Hspb1","Flt1","Esam","Eltd1","Abcb1a","Cldn5"), slot ="counts")
#Microglia
FeaturePlot(P28WT_Fixed.i, features= c("Csf1r" , "Cx3cr1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Csf1r" , "Cx3cr1"), slot ="counts")


#Selected markers Heatmap

DoHeatmap(P28WT_Fixed.i, features = c("Rasgrp","Smad3","Kcnip3","Mfsd4","Tanc1","Dsp","Rasgrf2","Gda","Rbm24"), slot ="scale.data")


#To be able to see the markers and not have a to wide range of expression we need to plot the log2(1+counts) expression
#(if we had 1 cell with an expression of 40 we were no able to see cells with expression of 15)

#Create a new assay to store log2(1+counts) information
#Save the linear normalized (LN) counts matrix
LN_counts_P28WT_Fixed.i<-as.data.frame(P28WT_Fixed.i@assays[["SCT"]]@counts)
#add +1 to each element of the matrix
LN_counts_P28WT_Fixed.i_1<-LN_counts_P28WT_Fixed.i+1
#check that it worked
LN_counts_P28WT_Fixed.i[1:10,1:4]
LN_counts_P28WT_Fixed.i_1[1:10,1:4]
#do the log2 from the counts+1 matrix
LN_counts_P28WT_Fixed.i_1_log2<-log2(LN_counts_P28WT_Fixed.i_1)
LN_counts_P28WT_Fixed.i_1_log2[1:10,1:4]
#Create a new assay to store log2(1+counts) information
log2_assay <- CreateAssayObject(counts = LN_counts_P28WT_Fixed.i_1_log2)

# add this assay to the previously created Seurat object
P28WT_Fixed.i[["LOG2"]] <- log2_assay

# Validate that the object now contains multiple assays
Assays(P28WT_Fixed.i)
#Check the active assay 
DefaultAssay(P28WT_Fixed.i)
# Switch the default to LOG2
DefaultAssay(P28WT_Fixed.i) <- "LOG2"
DefaultAssay(P28WT_Fixed.i)


#Granule mature 
FeaturePlot(P28WT_Fixed.i, features= c("Prox1","Plk5", "Ntng1", "Rbfox3", "Calb1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Prox1","Plk5", "Ntng1", "Rbfox3", "Calb1"), slot ="counts")
#Granule immature
FeaturePlot(P28WT_Fixed.i, features= c("Ncam1","Prox1","Dcx", "Fxyd7","Rbfox3", "Calb2" ), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Ncam1","Prox1","Dcx", "Fxyd7","Rbfox3", "Calb2" ), slot ="counts")
#Neuroblast 2
FeaturePlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Gal", "Sox11","Neurod1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Gal", "Sox11","Neurod1"), slot ="counts")
#Neuroblast1
FeaturePlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Eomes", "Igfbpl1", "Calb2", "Neurog2", "Sox11", "Neurod4", "Tac2", "Neurod1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Eomes", "Igfbpl1", "Calb2", "Neurog2", "Sox11", "Neurod4", "Tac2", "Neurod1"), slot ="counts")
#nIPC 
FeaturePlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Ascl", "Cdk1", "Tfap2c", "Eomes", "Igfbpl1","Neurog2", "Top2a", "Neurod4", "Aurkb", "Neurod1", "Ccnd2"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Prox1","Dcx","Ascl", "Cdk1", "Tfap2c", "Eomes", "Igfbpl1","Neurog2", "Top2a", "Neurod4", "Aurkb", "Neurod1", "Ccnd2"), slot ="counts")
#RGL
FeaturePlot(P28WT_Fixed.i, features= c("Tfap2c", "Sox9", "Lpar1", "Ascl1","Sox2", "Nes", "Vnn1","Hes5", "Rhcg", "Scl1a3", "Wnt8b", "Aldoc", "Apoe", "Id4", "Hopx", "GFAP"), slot ="counts")
VlnPlot(P28WT_Fixed.i, features= c("Tfap2c", "Sox9", "Lpar1", "Ascl1","Sox2", "Nes", "Vnn1","Hes5", "Rhcg", "Scl1a3", "Wnt8b", "Aldoc", "Apoe", "Id4", "Hopx", "GFAP"), slot ="counts")
#Astrocytes
FeaturePlot(P28WT_Fixed.i, features= c("S100b", "Fzd2", "Sox9","Sox2","Hes5","Scl1a3"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("S100b", "Fzd2", "Sox9","Sox2","Hes5","Scl1a3"), slot ="counts")
#OPC 
FeaturePlot(P28WT_Fixed.i, features= c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"), slot ="counts")
#Oligo 
FeaturePlot(P28WT_Fixed.i, features= c("Lpar1","Mbp", "Plp1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Lpar1","Mbp", "Plp1"), slot ="counts")
#Endothelial cells
FeaturePlot(P28WT_Fixed.i, features= c("Nes", "Prom1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Nes", "Prom1"), slot ="counts")
#Microglia
FeaturePlot(P28WT_Fixed.i, features= c("Csf1r" , "Cx3cr1"), slot ="counts",pt.size=1.5)
VlnPlot(P28WT_Fixed.i, features= c("Csf1r" , "Cx3cr1"), slot ="counts")

#------------------------------------------------Cluster investigation-----------------------------------------------------

#Playing with resolution
methanol.integrated <- FindNeighbors(methanol.integrated, dims = 1:12, k.param=10)
methanol.integrated<- FindClusters(methanol.integrated, resolution = seq(from = 0.1, to = 1, by = 0.1), algorithm=1)
clustree(methanol.integrated, prefix = "integrated_snn_res.")
#Playing with K
methanol.integrated_K5 <- FindNeighbors(methanol.integrated, dims = 1:12, k.param=5)
methanol.integrated_K7 <- FindNeighbors(methanol.integrated, dims = 1:12, k.param=7)
methanol.integrated_K10 <- FindNeighbors(methanol.integrated, dims = 1:12, k.param=10)
methanol.integrated_K15 <- FindNeighbors(methanol.integrated, dims = 1:12, k.param=15)
methanol.integrated_K20 <- FindNeighbors(methanol.integrated, dims = 1:12, k.param=20)

methanol.integrated_K5<- FindClusters(methanol.integrated_K5, resolution = 0.5, algorithm=1)
methanol.integrated_K7<- FindClusters(methanol.integrated_K7, resolution = 0.5, algorithm=1)
methanol.integrated_K10<- FindClusters(methanol.integrated_K10, resolution = 0.5, algorithm=1)
methanol.integrated_K15<- FindClusters(methanol.integrated_K15, resolution = 0.5, algorithm=1)
methanol.integrated_K20<- FindClusters(methanol.integrated_K20, resolution = 0.5, algorithm=1)



methanol.integrated[["K5"]] <- methanol.integrated_K5@meta.data[["seurat_clusters"]]
methanol.integrated[["K7"]] <- methanol.integrated_K7@meta.data[["seurat_clusters"]]
methanol.integrated[["K10"]] <- methanol.integrated_K10@meta.data[["seurat_clusters"]]
methanol.integrated[["K15"]] <- methanol.integrated_K15@meta.data[["seurat_clusters"]]
methanol.integrated[["K20"]] <- methanol.integrated_K20@meta.data[["seurat_clusters"]]

clustree(methanol.integrated, prefix = "K")

