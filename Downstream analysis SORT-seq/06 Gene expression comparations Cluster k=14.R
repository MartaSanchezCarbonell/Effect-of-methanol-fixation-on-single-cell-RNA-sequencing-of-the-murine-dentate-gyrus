#Setup the Seurat Object
#Install pakages (only once)
remove.packages("rlang")
install.packages("rlang")
install.packages('Seurat')
install.packages("devtools")
install.packages("umap-learn")
install.packages("xlsx")
install.packages("ggpubr")
install.packages("car")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dittoSeq")
BiocManager::install("edgeR")
install.packages("matrixStats")
install.packages("magrittr")
install.packages("Matrix")
install.packages("purrr")
install.packages("reshape2")
install.packages("S4Vectors")
install.packages("tibble")
install.packages("pheatmap")
install.packages("scales")
install.packages("ExperimentHub")
install.packages("pheatmap")
install.packages("Hmisc")
install.packages("MAST")
#Activate the following libraries/packages
library(rlang)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(devtools)
library(tidyverse)
library(xlsx)
library(dittoSeq)
library(matrixStats)
library(ggpubr)
library(car)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(pheatmap)
library(scales)
library(ExperimentHub)
library(pheatmap)
library(edgeR)
library(Hmisc)
library(MAST)
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

# chop of chromosome lables (Abel)
#chop_chr <- function  (name, splitcharacter = "__") {   strsplit(name, splitcharacter)[[1]][1]}
#pilot.data<-chop_chr(rownames(pilot.data), splitcharacter = "__")
#pilot.data


#when we create the seurat object we need to ask to remove the transcripts belonging to the ercc matrix and 
#add total_ercc_reads and percent.ERCC to object@meta.data in the total_ercc_reads and percent.ERCC column respectively
pilot <- CreateSeuratObject(counts = pilot.data[-ercc,], project = "Pilot analysis", min.cells = 0, min.features = 0)
pilot


pilot$orig.ident <- factor(pilot$orig.ident, levels=c("JEN-AU-s004","JEN-AU-s002"))

current.ids <- c("JEN-AU-s004","JEN-AU-s002")
new.ids <- c("Fresh_cells","Fixed_cells")
pilot@active.ident <- plyr::mapvalues(x = pilot@active.ident, from = current.ids, to = new.ids)
pilot$Condition <- Idents(pilot)

pilot$Condition <- factor(pilot$Condition, levels=c("Fresh_cells","Fixed_cells"))


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
pilot[["percent_ERCC"]] <- PercentageFeatureSet(pilot, pattern = "^ERCC")
pilot[["log10_nCount_RNA"]] <- log10(pilot$nCount_RNA)
pilot[["log10_nFeature_RNA"]] <- log10(pilot$nFeature_RNA)
pilot[["percent_ERCC_after_removal"]] <- PercentageFeatureSet(pilot, pattern = "^ERCC")
pilot[["ratio_mol_gen"]]<-pilot$nCount_RNA/pilot$nFeature_RNA
Control_of_ERCC_removal<- summary(pilot$percent_ERCC_after_removal)
Control_of_ERCC_removal


pilot.list<-SplitObject(pilot, split.by = "Condition")

Fresh_cells<-pilot.list$Fresh_cells
Fixed_cells<-pilot.list$Fixed_cells

#Histograms + cutoff lines

hist(pilot$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5), ylim=c(0,7))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fresh_cells$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5),ylim=c(0,7))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fixed_cells$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5),ylim=c(0,7))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")

hist(pilot$percent_mt, breaks=500, xlim = c(0,100))
hist(Fresh_cells$percent_mt,breaks=500, xlim = c(0,100))
hist(Fixed_cells$percent_mt,breaks=500, xlim = c(0,100))

hist(pilot$log10_nFeature_RNA, breaks=500, xlim = c(0,4))
abline(v=2.69897, col="blue")
hist(Fresh_cells$log10_nFeature_RNA,breaks=500, xlim = c(0,4))
abline(v=2.69897, col="blue")
hist(Fixed_cells$log10_nFeature_RNA,breaks=500, xlim = c(0,4))
abline(v=2.69897, col="blue")

hist(pilot$total_ercc_reads, breaks=500, xlim = c(0,2500))
abline(v=500, col="blue")
hist(Fresh_cells$total_ercc_reads,breaks=500, xlim = c(0,2500))
abline(v=500, col="blue")
hist(Fixed_cells$total_ercc_reads,breaks=500, xlim = c(0,2500))
abline(v=500, col="blue")

hist(pilot$ratio_mol_gen, breaks=500, xlim = c(0,55))
abline(v=1.2, col="blue")
hist(Fresh_cells$ratio_mol_gen,breaks=500, xlim = c(0,55))
abline(v=1.2, col="blue")
hist(Fixed_cells$ratio_mol_gen,breaks=500, xlim = c(0,55))
abline(v=1.2, col="blue")


#filtering
Fresh_cells<-subset(Fresh_cells, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
Fixed_cells<-subset(Fixed_cells, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
pilot<-subset(pilot, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
Fresh_cells
Fixed_cells
pilot


#Normalisation
Fresh_cells<-SCTransform(Fresh_cells,new.assay.name="SCT", variable.features.n=2000)
Fresh_cells

Fixed_cells<-SCTransform(Fixed_cells,new.assay.name="SCT", variable.features.n=2000)
Fixed_cells

#Integration 
#Create a list of Seurat objects 
pilot.list1<-list(Fresh_cells,Fixed_cells)
# select features for downstream integration
pilot.features <- SelectIntegrationFeatures(object.list = pilot.list1, nfeatures = 5000)
#run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated
pilot.list1<-PrepSCTIntegration(object.list = pilot.list1, anchor.features = pilot.features, 
                                verbose = FALSE)
#identify anchors
pilot.anchors<-FindIntegrationAnchors(object.list = pilot.list1, normalization.method = "SCT", 
                                      anchor.features = pilot.features, verbose = FALSE)
#integrate the datasets
pilot.integrated<-IntegrateData(anchorset = pilot.anchors, normalization.method = "SCT", 
                                verbose = FALSE, k.weight =64)

#once integrated we need to re-assign the order of the conditions
pilot.integrated$Condition <- factor(pilot.integrated$Condition, levels=c("Fresh_cells","Fixed_cells"))







#To analyze statistically significant differences between the gene levels of fresh and preserved cells : Number of significantly altered genes between fresh and fixed samples (at least 2-fol) significantly altered, with p-value < 0,05) (Wohnhaas 2019)







DefaultAssay(pilot.integrated) <- "RNA"
#Wilcoxon Rank Sum test

DE_pilot_wilcox<-FindMarkers(pilot.integrated,test.use="wilcox", slot="counts", logfc.threshold=0.3010299957, 
                             group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_pilot_wilcox

W_genes_altered<-nrow(DE_pilot_wilcox)  
W_genes_increased<-nrow(DE_pilot_wilcox[DE_pilot_wilcox$avg_log2FC > 0, ])
W_genes_decreased<-nrow(DE_pilot_wilcox[DE_pilot_wilcox$avg_log2FC < 0, ])


W_genes_altered
W_genes_increased
W_genes_decreased

DE_pilot_wilcox_summary<-data.frame(W_genes_altered, W_genes_increased, W_genes_decreased)
DE_pilot_wilcox_summary


#MAST

DE_pilot_MAST<-FindMarkers(pilot.integrated,test.use="MAST", slot="counts", logfc.threshold=0.3010299957, 
                             group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_pilot_MAST


M_genes_altered<-nrow(DE_pilot_MAST)  
M_genes_increased<-nrow(DE_pilot_MAST[DE_pilot_MAST$avg_log2FC > 0, ])
M_genes_decreased<-nrow(DE_pilot_MAST[DE_pilot_MAST$avg_log2FC < 0, ])


M_genes_altered
M_genes_increased
M_genes_decreased

DE_pilot_MAST_summary<-data.frame(M_genes_altered, M_genes_increased, M_genes_decreased)
DE_pilot_MAST_summary







#To gain further insights into preservation-related artifacts  Compared gene expression between preserved and freshly profiled samples in each cell type separately. (Denisenko 2020) 






DefaultAssay(pilot.integrated) <- "integrated"

#Perform linear dimensional reduction
pilot.integrated <- RunPCA(pilot.integrated, features = VariableFeatures(object = pilot.integrated))
# Examine and visualize PCA results a few different ways
DimPlot(pilot.integrated, reduction = "pca")
ElbowPlot(pilot.integrated)

#Cluster the cells (12 dims)
pilot.integrated <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=14)
pilot.integrated<- FindClusters(pilot.integrated, resolution = 0.5, algorithm=1)

#Lists of differential expressed genes (clusters)
pilot.markers <- FindAllMarkers(object=pilot.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox")
pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =2, order_by = avg_log2FC)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(pilot.markers, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                        paste(time.name), "pilot.markers 6 Gene expression comparations Cluster K=14.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)

#After cluster annotation with Linnarsson's database
#Assigning cell type identity to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("Astrocytes +  RGL", "nIPCs + Neuroblasts", "Immature neurons", 
                     "Mature neurons","Microglia","Oligodendrocytes", "OPCs", "Endothelial cells")
pilot.integrated@active.ident <- plyr::mapvalues(x = pilot.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
pilot.integrated[["new.cluster.ids"]]<-pilot.integrated@active.ident

DefaultAssay(pilot.integrated) <- "RNA"

#Subset cell identities
Astrocytes_RGL<-subset(pilot.integrated, subset = new.cluster.ids=="Astrocytes +  RGL")
nIPCs_Neuroblasts<-subset(pilot.integrated, subset = new.cluster.ids=="nIPCs + Neuroblasts")
Immature_neurons<-subset(pilot.integrated, subset = new.cluster.ids=="Immature neurons")
Mature_neurons<-subset(pilot.integrated, subset = new.cluster.ids=="Mature neurons")
Microglia<-subset(pilot.integrated, subset = new.cluster.ids=="Microglia")
Oligodendrocytes<-subset(pilot.integrated, subset = new.cluster.ids=="Oligodendrocytes")
OPCs<-subset(pilot.integrated, subset = new.cluster.ids=="OPCs")
Endothelial_cells<-subset(pilot.integrated, subset = new.cluster.ids=="Endothelial cells")

#Differential gene expression tests

#Wilcoxon Rank Sum test
#Astrocytes_RGL
DE_Astrocytes_RGL_wilcox<-FindMarkers(Astrocytes_RGL,test.use="wilcox", slot="counts", min.pct = 0.5,min.diff.pct=0.05, logfc.threshold=1, 
                             group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_Astrocytes_RGL_wilcox

W_Astrocytes_RGL_altered<-nrow(DE_Astrocytes_RGL_wilcox)  
W_Astrocytes_RGL_genes_increased<-nrow(DE_Astrocytes_RGL_wilcox[DE_Astrocytes_RGL_wilcox$avg_log2FC > 0, ])
W_Astrocytes_RGL_genes_decreased<-nrow(DE_Astrocytes_RGL_wilcox[DE_Astrocytes_RGL_wilcox$avg_log2FC < 0, ])

DE_Astrocytes_RGL_wilcox_summary<-data.frame(W_Astrocytes_RGL_altered, W_Astrocytes_RGL_genes_increased, W_Astrocytes_RGL_genes_decreased)
DE_Astrocytes_RGL_wilcox_summary


#nIPCs_Neuroblasts
DE_nIPCs_Neuroblasts_wilcox<-FindMarkers(nIPCs_Neuroblasts,test.use="wilcox", slot="counts", min.pct = 0.5,min.diff.pct=0.05, logfc.threshold=1, 
                                      group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_nIPCs_Neuroblasts_wilcox

W_nIPCs_Neuroblasts_altered<-nrow(DE_nIPCs_Neuroblasts_wilcox)  
W_nIPCs_Neuroblasts_genes_increased<-nrow(DE_nIPCs_Neuroblasts_wilcox[DE_nIPCs_Neuroblasts_wilcox$avg_log2FC > 0, ])
W_nIPCs_Neuroblasts_genes_decreased<-nrow(DE_nIPCs_Neuroblasts_wilcox[DE_nIPCs_Neuroblasts_wilcox$avg_log2FC < 0, ])

DE_nIPCs_Neuroblasts_wilcox_summary<-data.frame(W_nIPCs_Neuroblasts_altered, W_nIPCs_Neuroblasts_genes_increased, W_nIPCs_Neuroblasts_genes_decreased)
DE_nIPCs_Neuroblasts_wilcox_summary


#Immature_neurons
DE_Immature_neurons_wilcox<-FindMarkers(Immature_neurons,test.use="wilcox", slot="counts", min.pct = 0.5,min.diff.pct=0.05, logfc.threshold=1, 
                                         group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_Immature_neurons_wilcox

W_Immature_neurons_altered<-nrow(DE_Immature_neurons_wilcox)  
W_Immature_neurons_genes_increased<-nrow(DE_Immature_neurons_wilcox[DE_Immature_neurons_wilcox$avg_log2FC > 0, ])
W_Immature_neurons_genes_decreased<-nrow(DE_Immature_neurons_wilcox[DE_Immature_neurons_wilcox$avg_log2FC < 0, ])

DE_Immature_neurons_wilcox_summary<-data.frame(W_Immature_neurons_altered, W_Immature_neurons_genes_increased, W_Immature_neurons_genes_decreased)
DE_Immature_neurons_wilcox_summary

#Mature_neurons
DE_Mature_neurons_wilcox<-FindMarkers(Mature_neurons,test.use="wilcox", slot="counts", min.pct = 0.5,min.diff.pct=0.05, logfc.threshold=1, 
                                        group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_Mature_neurons_wilcox

W_Mature_neurons_altered<-nrow(DE_Mature_neurons_wilcox)  
W_Mature_neurons_genes_increased<-nrow(DE_Mature_neurons_wilcox[DE_Mature_neurons_wilcox$avg_log2FC > 0, ])
W_Mature_neurons_genes_decreased<-nrow(DE_Mature_neurons_wilcox[DE_Mature_neurons_wilcox$avg_log2FC < 0, ])

DE_Mature_neurons_wilcox_summary<-data.frame(W_Mature_neurons_altered, W_Mature_neurons_genes_increased, W_Mature_neurons_genes_decreased)
DE_Mature_neurons_wilcox_summary

