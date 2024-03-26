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
install.packages('parallel')
install.packages("ROCR")
install.packages("KernSmooth")
install.packages("fields")
install.packages("Matrix")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
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
library(parallel)
library(ROCR)
library(KernSmooth)
library(fields)
library(Matrix)
library(DoubletFinder)
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

#pilot[["Stmn2"]]<-FetchData(object = pilot, vars = "Stmn2")
#pilot[["Mog"]]<-FetchData(object = pilot, vars = "Mog")
#pilot[["Aldoc"]]<-FetchData(object = pilot, vars = "Aldoc")
#pilot[["C1qc"]]<-FetchData(object = pilot, vars = "C1qc")
#pilot[["Cldn5"]]<-FetchData(object = pilot, vars = "Cldn5")



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
Fresh_cells_filt<-subset(Fresh_cells, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
Fixed_cells_filt<-subset(Fixed_cells, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
pilot_filt<-subset(pilot, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
pilot_filt
Fresh_cells_filt
Fixed_cells_filt




              #Hochgerner 2018 





#remove doublet cells based on coexpression (>1 molecule) of any pair of the following marker genes: Stmn2 (neurons), Mog (oligodendrocytes), Aldoc (astrocytes), C1qc (microglia), and Cldn5 (endothelial).

#Stmn2<- grep("Stmn2",pilot_filt@assays[["RNA"]]@data@Dimnames[[1]], value=TRUE) 
#Stmn2
#Mog<- grep("Mog",pilot_filt@assays[["RNA"]]@data@Dimnames[[1]], value=TRUE)
#Mog
#Aldoc<- grep("Aldoc",pilot_filt@assays[["RNA"]]@data@Dimnames[[1]], value=TRUE)
#Aldoc
#C1qc<- grep("C1qc",pilot_filt@assays[["RNA"]]@data@Dimnames[[1]], value=TRUE)
#C1qc
#Cldn5<- grep("Cldn5",pilot_filt@assays[["RNA"]]@data@Dimnames[[1]], value=TRUE)
#Cldn5
 
#Stmn2 <- FetchData(object = pilot_filt, vars = "Stmn2")
#Mog<- FetchData(object = pilot_filt, vars = "Mog")
#Aldoc<- FetchData(object = pilot_filt, vars = "Aldoc")
#C1qc<- FetchData(object = pilot_filt, vars = "C1qc")
#Cldn5<- FetchData(object = pilot_filt, vars = "Cldn5")

#With this we subset the cells that are fulfilling this criteria

H_Doublets_pilot_filt<-subset(pilot_filt, subset = rna_Stmn2 > 1 & rna_Mog>1 | rna_Stmn2>1 & rna_Aldoc>1 |rna_Stmn2>1 & rna_C1qc>1 | 
                              rna_Stmn2>1 & rna_Cldn5>1 |rna_Mog>1 & rna_Aldoc>1 |rna_Mog>1 & rna_C1qc>1 |rna_Mog>1 & rna_Cldn5>1 |
                              rna_Aldoc>1 & rna_C1qc>1 |rna_Aldoc>1 & rna_Cldn5>1 | rna_C1qc>1 & rna_Cldn5>1)

H_Doublets_Fresh_cells_filt<-subset(Fresh_cells_filt, subset = rna_Stmn2 > 1 & rna_Mog>1 | rna_Stmn2>1 & rna_Aldoc>1 |rna_Stmn2>1 & rna_C1qc>1 | 
                                    rna_Stmn2>1 & rna_Cldn5>1 |rna_Mog>1 & rna_Aldoc>1 |rna_Mog>1 & rna_C1qc>1 |rna_Mog>1 & rna_Cldn5>1 |
                                    rna_Aldoc>1 & rna_C1qc>1 |rna_Aldoc>1 & rna_Cldn5>1 | rna_C1qc>1 & rna_Cldn5>1)

H_Doublets_Fixed_cells_filt<-subset(Fixed_cells_filt, subset = rna_Stmn2 > 1 & rna_Mog>1 | rna_Stmn2>1 & rna_Aldoc>1 |rna_Stmn2>1 & rna_C1qc>1 | 
                                    rna_Stmn2>1 & rna_Cldn5>1 |rna_Mog>1 & rna_Aldoc>1 |rna_Mog>1 & rna_C1qc>1 |rna_Mog>1 & rna_Cldn5>1 |
                                    rna_Aldoc>1 & rna_C1qc>1 |rna_Aldoc>1 & rna_Cldn5>1 | rna_C1qc>1 & rna_Cldn5>1)



#H_Doublets_pilot_filt<-pilot_filt[, which(x = Stmn2 > 1 & Mog>1 | Stmn2>1 & Aldoc>1 |Stmn2>1 & C1qc>1 | 
               #Stmn2>1 & Cldn5>1 |Mog>1 & Aldoc>1 |Mog>1 & C1qc>1 |Mog>1 & Cldn5>1 |
                #Aldoc>1 & C1qc>1 |Aldoc>1 & Cldn5>1 | C1qc>1 & Cldn5>1)]


#H_Doublets_Fresh_cells_filt<-Fresh_cells_filt[, which(x = Stmn2 > 1 & Mog>1 | Stmn2>1 & Aldoc>1 |Stmn2>1 & C1qc>1 | 
                                                      #Stmn2>1 & Cldn5>1 |Mog>1 & Aldoc>1 |Mog>1 & C1qc>1 |Mog>1 & Cldn5>1 | 
                                                      #Aldoc>1 & C1qc>1 |Aldoc>1 & Cldn5>1 | C1qc>1 & Cldn5>1)]

#H_Doublets_Fixed_cells_filt<-Fixed_cells_filt[, which(x = Stmn2 > 1 & Mog>1 | Stmn2>1 & Aldoc>1 |Stmn2>1 & C1qc>1 | 
                                          #Stmn2>1 & Cldn5>1 |Mog>1 & Aldoc>1 |Mog>1 & C1qc>1 |Mog>1 & Cldn5>1 |
                                          #Aldoc>1 & C1qc>1 |Aldoc>1 & Cldn5>1 | C1qc>1 & Cldn5>1)]


pilot_filt
Fresh_cells_filt
Fixed_cells_filt

H_Doublets_pilot_filt
H_Doublets_Fresh_cells_filt
H_Doublets_Fixed_cells_filt

colnames(H_Doublets_pilot_filt)
colnames(H_Doublets_Fresh_cells_filt)
colnames(H_Doublets_Fixed_cells_filt)

H_Doublets_pilot_filt@assays[["RNA"]][c("Stmn2","Mog","Aldoc","C1qc","Cldn5"),]
H_Doublets_Fresh_cells_filt@assays[["RNA"]][c("Stmn2","Mog","Aldoc","C1qc","Cldn5"),]
H_Doublets_Fixed_cells_filt@assays[["RNA"]][c("Stmn2","Mog","Aldoc","C1qc","Cldn5"),]

H_Doublets_pilot_filt@assays[["RNA"]]@scale.data[c("Stmn2","Mog","Aldoc","C1qc","Cldn5"),]
H_Doublets_Fresh_cells_filt@assays[["RNA"]]@scale.data[c("Stmn2","Mog","Aldoc","C1qc","Cldn5"),]
H_Doublets_Fixed_cells_filt@assays[["RNA"]]@scale.data[c("Stmn2","Mog","Aldoc","C1qc","Cldn5"),]



#Histograms of the doublets
      #Pilot
hist(H_Doublets_pilot_filt$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5), ylim=c(0,7))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")

hist(H_Doublets_pilot_filt$percent_mt, breaks=500, xlim = c(0,100))

hist(H_Doublets_pilot_filt$log10_nFeature_RNA, breaks=500, xlim = c(0,4))
abline(v=2.69897, col="blue")

hist(H_Doublets_pilot_filt$total_ercc_reads, breaks=500, xlim = c(0,2500))
abline(v=500, col="blue")

hist(H_Doublets_pilot_filt$ratio_mol_gen, breaks=500, xlim = c(0,55))
abline(v=1.2, col="blue")

      #Fresh
hist(H_Doublets_Fresh_cells_filt$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5), ylim=c(0,7))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")

hist(H_Doublets_Fresh_cells_filt$percent_mt, breaks=500, xlim = c(0,100))

hist(H_Doublets_Fresh_cells_filt$log10_nFeature_RNA, breaks=500, xlim = c(0,4))
abline(v=2.69897, col="blue")

hist(H_Doublets_Fresh_cells_filt$total_ercc_reads, breaks=500, xlim = c(0,2500))
abline(v=500, col="blue")

hist(H_Doublets_Fresh_cells_filt$ratio_mol_gen, breaks=500, xlim = c(0,55))
abline(v=1.2, col="blue")

      #Fixed
hist(H_Doublets_Fixed_cells_filt$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5), ylim=c(0,7))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")

hist(H_Doublets_Fixed_cells_filt$percent_mt, breaks=500, xlim = c(0,100))

hist(H_Doublets_Fixed_cells_filt$log10_nFeature_RNA, breaks=500, xlim = c(0,4))
abline(v=2.69897, col="blue")

hist(H_Doublets_Fixed_cells_filt$total_ercc_reads, breaks=500, xlim = c(0,2500))
abline(v=500, col="blue")

hist(H_Doublets_Fixed_cells_filt$ratio_mol_gen, breaks=500, xlim = c(0,55))
abline(v=1.2, col="blue")



#Normalising
H_Doublets_pilot_filt <- NormalizeData(H_Doublets_pilot_filt, normalization.method = "LogNormalize", scale.factor= TRUE, margin = TRUE)
H_Doublets_Fresh_cells_filt <- NormalizeData(H_Doublets_Fresh_cells_filt, normalization.method = "LogNormalize", scale.factor= TRUE, margin = TRUE)
H_Doublets_Fixed_cells_filt <- NormalizeData(H_Doublets_Fixed_cells_filt, normalization.method = "LogNormalize", scale.factor= TRUE, margin = TRUE)

#Identification of highly variable features (feature selection)
H_Doublets_pilot_filt <- FindVariableFeatures(H_Doublets_pilot_filt, selection.method = "vst", nfeatures = 2000)
H_Doublets_Fresh_cells_filt <- FindVariableFeatures(H_Doublets_Fresh_cells_filt, selection.method = "vst", nfeatures = 2000)
H_Doublets_Fixed_cells_filt <- FindVariableFeatures(H_Doublets_Fixed_cells_filt, selection.method = "vst", nfeatures = 2000)


#Scaling the data
all.genes_H_Doublets_pilot_filt <- rownames(H_Doublets_pilot_filt)
H_Doublets_pilot_filt <- ScaleData(H_Doublets_pilot_filt, features = all.genes_H_Doublets_pilot_filt)

all.genes_H_Doublets_Fresh_cells_filt <- rownames(H_Doublets_Fresh_cells_filt)
H_Doublets_Fresh_cells_filt <- ScaleData(H_Doublets_Fresh_cells_filt, features = all.genes_H_Doublets_Fresh_cells_filt)

all.genes_H_Doublets_Fixed_cells_filt <- rownames(H_Doublets_Fixed_cells_filt)
H_Doublets_Fixed_cells_filt <- ScaleData(H_Doublets_Fixed_cells_filt, features = all.genes_H_Doublets_Fixed_cells_filt)

#Perform linear dimensional reduction
H_Doublets_pilot_filt <- RunPCA(H_Doublets_pilot_filt, features = VariableFeatures(object = H_Doublets_pilot_filt), npcs = 10)
H_Doublets_Fresh_cells_filt <- RunPCA(H_Doublets_Fresh_cells_filt, features = VariableFeatures(object = H_Doublets_Fresh_cells_filt),npcs = 2)
H_Doublets_Fixed_cells_filt <- RunPCA(H_Doublets_Fixed_cells_filt, features = VariableFeatures(object = H_Doublets_Fixed_cells_filt),npcs = 7)
# Examine and visualize PCA results a few different ways
#DimPlot(Doublets, reduction = "pca")
#Determine the ‘dimensionality’ of the dataset
#Doublets <- JackStraw(Doublets, num.replicate = 100)
#Doublets <- ScoreJackStraw(Doublets, dims = 1:20)
#JackStrawPlot(Doublets, dims = 1:20)
ElbowPlot(H_Doublets_pilot_filt)
ElbowPlot(H_Doublets_Fresh_cells_filt)
ElbowPlot(H_Doublets_Fixed_cells_filt)


#Cluster the cells (18 dims)
H_Doublets_pilot_filt <- FindNeighbors(H_Doublets_pilot_filt, dims = 1:10, k.param=10)
H_Doublets_pilot_filt <- FindClusters(H_Doublets_pilot_filt, resolution = 0.5, algorithm=1)

H_Doublets_Fresh_cells_filt <- FindNeighbors(H_Doublets_Fresh_cells_filt, dims = 1:2, k.param=10)
H_Doublets_Fresh_cells_filt <- FindClusters(H_Doublets_Fresh_cells_filt, resolution = 0.5, algorithm=1)

H_Doublets_Fixed_cells_filt <- FindNeighbors(H_Doublets_Fixed_cells_filt, dims = 1:7, k.param=10)
H_Doublets_Fixed_cells_filt <- FindClusters(H_Doublets_Fixed_cells_filt, resolution = 0.5, algorithm=1)
#Run non-linear dimensional reduction (UMAP/tSNE)


#UMAP
H_Doublets_pilot_filt <- RunUMAP(H_Doublets_pilot_filt, dims = 1:10, n.neighbors=10L)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(H_Doublets_pilot_filt, reduction = "umap")

H_Doublets_Fresh_cells_filt <- RunUMAP(H_Doublets_Fresh_cells_filt, dims = 1:2, n.neighbors=2L)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(H_Doublets_Fresh_cells_filt, reduction = "umap")

H_Doublets_Fixed_cells_filt <- RunUMAP(H_Doublets_Fixed_cells_filt, dims = 1:7, n.neighbors=7L)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(H_Doublets_Fixed_cells_filt, reduction = "umap")


# Counts

FeaturePlot(H_Doublets_pilot_filt, features = c("Stmn2", "Mog"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)
FeaturePlot(H_Doublets_pilot_filt, features = c("Stmn2", "Aldoc"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("Stmn2", "C1qc"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("Stmn2", "Cldn5"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)
FeaturePlot(H_Doublets_pilot_filt, features = c("Mog", "Aldoc"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("Mog", "C1qc"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("Mog", "Cldn5"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)
FeaturePlot(H_Doublets_pilot_filt, features = c("Aldoc", "C1qc"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)
FeaturePlot(H_Doublets_pilot_filt, features = c("Aldoc", "Cldn5"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("C1qc", "Cldn5"), blend = TRUE, pt.size=7, slot="counts",blend.threshold=0)

#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Stmn2", "Mog"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Stmn2", "Aldoc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Stmn2", "C1qc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Stmn2", "Cldn5"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Mog", "Aldoc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Mog", "C1qc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Mog", "Cldn5"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Aldoc", "C1qc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Aldoc", "Cldn5"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("C1qc", "Cldn5"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)

FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Stmn2", "Mog"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Stmn2", "Aldoc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Stmn2", "C1qc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Stmn2", "Cldn5"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Mog", "Aldoc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Mog", "C1qc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Mog", "Cldn5"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Aldoc", "C1qc"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Aldoc", "Cldn5"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("C1qc", "Cldn5"), blend = TRUE, pt.size=7, slot="counts", blend.threshold=0)


#scaled.data
FeaturePlot(H_Doublets_pilot_filt, features = c("Stmn2", "Mog"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)
FeaturePlot(H_Doublets_pilot_filt, features = c("Stmn2", "Aldoc"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("Stmn2", "C1qc"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("Stmn2", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)
FeaturePlot(H_Doublets_pilot_filt, features = c("Mog", "Aldoc"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("Mog", "C1qc"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("Mog", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)
FeaturePlot(H_Doublets_pilot_filt, features = c("Aldoc", "C1qc"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)
FeaturePlot(H_Doublets_pilot_filt, features = c("Aldoc", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)
#FeaturePlot(H_Doublets_pilot_filt, features = c("C1qc", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data",blend.threshold=0)

#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Stmn2", "Mog"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Stmn2", "Aldoc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Stmn2", "C1qc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Stmn2", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Mog", "Aldoc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Mog", "C1qc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Mog", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Aldoc", "C1qc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("Aldoc", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fresh_cells_filt, features = c("C1qc", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)

FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Stmn2", "Mog"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Stmn2", "Aldoc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Stmn2", "C1qc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Stmn2", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Mog", "Aldoc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Mog", "C1qc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Mog", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Aldoc", "C1qc"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("Aldoc", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)
#FeaturePlot(H_Doublets_Fixed_cells_filt, features = c("C1qc", "Cldn5"), blend = TRUE, pt.size=7, slot="scale.data", blend.threshold=0)











#Predict doublets with Seurat







#Here, we will use DoubletFinder to predict doublet cells. But before doing doublet detection we need to run scaling, 
#variable gene selection and pca, as well as UMAP for visualization. 
pilot_filt
Fresh_cells_filt
Fixed_cells_filt


#Normalisation
pilot_filt<-SCTransform(pilot_filt,new.assay.name="SCT", variable.features.n=2000)
pilot_filt

Fresh_cells_filt<-SCTransform(Fresh_cells_filt,new.assay.name="SCT", variable.features.n=2000)
Fresh_cells_filt

Fixed_cells_filt<-SCTransform(Fixed_cells_filt,new.assay.name="SCT", variable.features.n=2000)
Fixed_cells_filt


#Integration 
#Create a list of Seurat objects 
pilot.list1<-list(Fresh_cells_filt,Fixed_cells_filt)
# select features for downstream integration
pilot.features <- SelectIntegrationFeatures(object.list = pilot.list1, nfeatures = 2000)
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
pilot.integrated$Condition <- factor(pilot.integrated$Condition, levels=c("Fresh_cells_filt","Fixed_cells_filt"))


#Perform linear dimensional reduction
pilot_filt <- RunPCA(pilot_filt, features = VariableFeatures(object = pilot_filt))
pilot.integrated <- RunPCA(pilot.integrated, features = VariableFeatures(object = pilot_filt))
#Fresh_cells_filt <- RunPCA(Fresh_cells_filt, features = VariableFeatures(object = Fresh_cells_filt))
#Fixed_cells_filt <- RunPCA(Fixed_cells_filt, features = VariableFeatures(object = Fixed_cells_filt))

ElbowPlot(pilot_filt)
ElbowPlot(pilot.integrated)
#ElbowPlot(Fresh_cells_filt)
#ElbowPlot(Fixed_cells_filt)


#Cluster the cells (18 dims)
pilot_filt <- FindNeighbors(pilot_filt, dims = 1:12, k.param=10)
pilot_filt <- FindClusters(pilot_filt, resolution = 0.5, algorithm=1)

pilot.integrated <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=10)
pilot.integrated <- FindClusters(pilot.integrated, resolution = 0.5, algorithm=1)

#Fresh_cells_filt <- FindNeighbors(Fresh_cells_filt, dims = 1:12, k.param=10)
#Fresh_cells_filt <- FindClusters(Fresh_cells_filt, resolution = 0.5, algorithm=1)

#Fixed_cells_filt <- FindNeighbors(Fixed_cells_filt, dims = 1:12, k.param=10)
#Fixed_cells_filt <- FindClusters(Fixed_cells_filt, resolution = 0.5, algorithm=1)
#Run non-linear dimensional reduction (UMAP/tSNE)


#UMAP
pilot_filt <- RunUMAP(pilot_filt, dims = 1:12)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pilot_filt, reduction = "umap")

pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:12)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pilot.integrated, reduction = "umap")

#Fresh_cells_filt <- RunUMAP(Fresh_cells_filt, dims = 1:12)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
#DimPlot(Fresh_cells_filt, reduction = "umap")

#Fixed_cells_filt <- RunUMAP(Fixed_cells_filt, dims = 1:12)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
#DimPlot(Fixed_cells_filt, reduction = "umap")


              #pilot_filt

#DoubletFinder
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_pilot_filt <- paramSweep_v3(pilot_filt, PCs = 1:10, sct = TRUE)
sweep.stats_pilot_filt <- summarizeSweep(sweep.res.list_pilot_filt, GT = FALSE)
bcmvn_pilot_filt <- find.pK(sweep.stats_pilot_filt)
bcmvn_pilot_filt

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(pilot_filt$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.004*nrow(pilot@meta.data))  ## Assuming 0.4% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
homotypic.prop
nExp_poi
nExp_poi.adj


## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
pilot_filt <- doubletFinder_v3(pilot_filt, PCs = 1:10, pN = 0.25, pK = 0.16, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
pilot_filt <- doubletFinder_v3(pilot_filt, PCs = 1:10, pN = 0.25, pK = 0.16, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.16_3", sct = TRUE)
DimPlot(pilot_filt, reduction = "umap",group.by="DF.classifications_0.25_0.16_3", pt.size=1.5)


#To see how many cells are doublets and especifically which ones
Doublets_pilot<-subset(pilot_filt, DF.classifications_0.25_0.16_3 == "Doublet")
  
Doublets_pilot
colnames(Doublets_pilot)



              #pilot.integrated

#DoubletFinder
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_pilot.integrated <- paramSweep_v3(pilot.integrated, PCs = 1:10, sct = TRUE)
sweep.stats_pilot.integrated <- summarizeSweep(sweep.res.list_pilot.integrated, GT = FALSE)
bcmvn_pilot.integrated <- find.pK(sweep.stats_pilot.integrated)
bcmvn_pilot.integrated

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop_pilot.integrated <- modelHomotypic(pilot.integrated$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi_pilot.integrated <- round(0.004*nrow(pilot@meta.data))  ## Assuming 0.4% doublet formation rate - tailor for your dataset
nExp_poi.adj_pilot.integrated <- round(nExp_poi_pilot.integrated*(1-homotypic.prop_pilot.integrated))
homotypic.prop_pilot.integrated
nExp_poi_pilot.integrated
nExp_poi.adj_pilot.integrated


## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
pilot.integrated <- doubletFinder_v3(pilot.integrated, PCs = 1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi_pilot.integrated, reuse.pANN = FALSE, sct = TRUE)
pilot.integrated <- doubletFinder_v3(pilot.integrated, PCs = 1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi.adj_pilot.integrated, reuse.pANN = "pANN_0.25_0.06_3", sct = TRUE)
DimPlot(pilot.integrated, reduction = "umap",group.by="DF.classifications_0.25_0.06_3", pt.size=1.5)


#To see how many cells are doublets and especifically which ones
Doublets_pilot.integrated<-subset(pilot.integrated, DF.classifications_0.25_0.06_3 == "Doublet")

Doublets_pilot.integrated
colnames(Doublets_pilot.integrated)











#DoubletDecon





#devtools::install_github('EDePasquale/DoubletDecon')



