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
BiocManager::install("limma")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
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
install.packages('enrichR')
install.packages("ggpubr")
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
library(limma)
library(enrichR)
library(ggpubr)
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
                             group.by="Condition", ident.1="Fresh_cells", ident.2="Fixed_cells")
DE_pilot_wilcox

W_genes_altered<-nrow(DE_pilot_wilcox)  
W_genes_increased<-nrow(DE_pilot_wilcox[DE_pilot_wilcox$avg_log2FC > 0, ])
W_genes_decreased<-nrow(DE_pilot_wilcox[DE_pilot_wilcox$avg_log2FC < 0, ])


W_genes_altered
W_genes_increased
W_genes_decreased

DE_pilot_wilcox_summary<-data.frame(W_genes_altered, W_genes_increased, W_genes_decreased)
DE_pilot_wilcox_summary

#Saveing the marker list to run GO anlysis
time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(DE_pilot_wilcox, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "DE_Analysis Fresh vs Fixed.xlsx"))
           , sheetName = "DE Genes", col.names = TRUE, row.names = TRUE, append = FALSE)



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
pilot.integrated <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=10)
pilot.integrated<- FindClusters(pilot.integrated, resolution = 0.5, algorithm=1)

#Lists of differential expressed genes (clusters)
pilot.markers <- FindAllMarkers(object=pilot.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox")
pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =2, order_by = avg_log2FC)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(pilot.markers, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                        paste(time.name), "pilot.markers 6 Gene expression comparations.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)

#After cluster annotation with Linnarsson's database
#Assigning cell type identity to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("Astrocytes +  RGL", "nIPCs + Neuroblasts", "Neurons", 
                     "Neurons","Microglia","OPCs + Oligodendrocytes", "OPCs + Oligodendrocytes", "Endothelial cells")
pilot.integrated@active.ident <- plyr::mapvalues(x = pilot.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
pilot.integrated[["new.cluster.ids"]]<-pilot.integrated@active.ident

DefaultAssay(pilot.integrated) <- "RNA"

#Subset cell identities
Astrocytes_RGL<-subset(pilot.integrated, subset = new.cluster.ids=="Astrocytes +  RGL")
nIPCs_Neuroblasts<-subset(pilot.integrated, subset = new.cluster.ids=="nIPCs + Neuroblasts")
Neurons<-subset(pilot.integrated, subset = new.cluster.ids=="Neurons")
Microglia<-subset(pilot.integrated, subset = new.cluster.ids=="Microglia")
OPCs_Oligodendrocytes<-subset(pilot.integrated, subset = new.cluster.ids=="OPCs + Oligodendrocytes")
Endothelial_cells<-subset(pilot.integrated, subset = new.cluster.ids=="Endothelial cells")

#number of cells per cell identity
table(Idents(pilot.integrated), pilot.integrated$Condition)






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


#Neurons
DE_Neurons_wilcox<-FindMarkers(Neurons,test.use="wilcox", slot="counts", min.pct = 0.5,min.diff.pct=0.05, logfc.threshold=1, 
                                         group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_Neurons_wilcox

W_Neurons_altered<-nrow(DE_Neurons_wilcox)  
W_Neurons_genes_increased<-nrow(DE_Neurons_wilcox[DE_Neurons_wilcox$avg_log2FC > 0, ])
W_Neurons_genes_decreased<-nrow(DE_Neurons_wilcox[DE_Neurons_wilcox$avg_log2FC < 0, ])

DE_Neurons_wilcox_summary<-data.frame(W_Neurons_altered, W_Neurons_genes_increased, W_Neurons_genes_decreased)
DE_Neurons_wilcox_summary

#Microglia
DE_Microglia_wilcox<-FindMarkers(Microglia,test.use="wilcox", slot="counts", min.pct = 0.5,min.diff.pct=0.05, logfc.threshold=1, 
                                        group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_Microglia_wilcox

W_Microglia_altered<-nrow(DE_Microglia_wilcox)  
W_Microglia_genes_increased<-nrow(DE_Microglia_wilcox[DE_Microglia_wilcox$avg_log2FC > 0, ])
W_Microglia_genes_decreased<-nrow(DE_Microglia_wilcox[DE_Microglia_wilcox$avg_log2FC < 0, ])

DE_Microglia_wilcox_summary<-data.frame(W_Microglia_altered, W_Microglia_genes_increased, W_Microglia_genes_decreased)
DE_Microglia_wilcox_summary

#OPCs_Oligodendrocytes
DE_OPCs_Oligodendrocytes_wilcox<-FindMarkers(OPCs_Oligodendrocytes,test.use="wilcox", slot="counts", min.pct = 0.5,min.diff.pct=0.05, logfc.threshold=1, 
                                 group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_OPCs_Oligodendrocytes_wilcox

W_OPCs_Oligodendrocytes_altered<-nrow(DE_OPCs_Oligodendrocytes_wilcox)  
W_OPCs_Oligodendrocytes_genes_increased<-nrow(DE_OPCs_Oligodendrocytes_wilcox[DE_OPCs_Oligodendrocytes_wilcox$avg_log2FC > 0, ])
W_OPCs_Oligodendrocytes_genes_decreased<-nrow(DE_OPCs_Oligodendrocytes_wilcox[DE_OPCs_Oligodendrocytes_wilcox$avg_log2FC < 0, ])

DE_OPCs_Oligodendrocytes_wilcox_summary<-data.frame(W_OPCs_Oligodendrocytes_altered, W_OPCs_Oligodendrocytes_genes_increased, W_OPCs_Oligodendrocytes_genes_decreased)
DE_OPCs_Oligodendrocytes_wilcox_summary

#Endothelial_cells
DE_Endothelial_cells_wilcox<-FindMarkers(Endothelial_cells,test.use="wilcox", slot="counts", min.pct = 0.5,min.diff.pct=0.05, logfc.threshold=1, 
                                             group.by="Condition", ident.1="Fixed_cells", ident.2="Fresh_cells")
DE_Endothelial_cells_wilcox

W_Endothelial_cells_altered<-nrow(DE_Endothelial_cells_wilcox)  
W_Endothelial_cells_genes_increased<-nrow(DE_Endothelial_cells_wilcox[DE_Endothelial_cells_wilcox$avg_log2FC > 0, ])
W_Endothelial_cells_genes_decreased<-nrow(DE_Endothelial_cells_wilcox[DE_Endothelial_cells_wilcox$avg_log2FC < 0, ])

DE_Endothelial_cells_wilcox_summary<-data.frame(W_Endothelial_cells_altered, W_Endothelial_cells_genes_increased, W_Endothelial_cells_genes_decreased)
DE_Endothelial_cells_wilcox_summary

#Summary
Indentities<-c("Astrocytes +  RGL", "nIPCs + Neuroblasts", 
               "Neurons","Microglia", "OPCs + Oligodendrocytes", "Endothelial cells")
Altered_genes<-c(W_Astrocytes_RGL_altered,W_nIPCs_Neuroblasts_altered,W_Neurons_altered,
                 W_Microglia_altered,W_OPCs_Oligodendrocytes_altered,W_Endothelial_cells_altered)
Increased_genes<-c(W_Astrocytes_RGL_genes_increased,W_nIPCs_Neuroblasts_genes_increased,W_Neurons_genes_increased,
                 W_Microglia_genes_increased,W_OPCs_Oligodendrocytes_genes_increased,W_Endothelial_cells_genes_increased)
Decreased_genes<-c(W_Astrocytes_RGL_genes_decreased,W_nIPCs_Neuroblasts_genes_decreased,W_Neurons_genes_decreased,
                   W_Microglia_genes_decreased,W_OPCs_Oligodendrocytes_genes_decreased,W_Endothelial_cells_genes_decreased)
W_summary<-data.frame(Indentities,Altered_genes,Increased_genes,Decreased_genes)
W_summary

#GO Analysis

Idents(Astrocytes_RGL) <- "Condition"
Idents(nIPCs_Neuroblasts) <- "Condition"
Idents(Neurons) <- "Condition"
Idents(Microglia) <- "Condition"
Idents(OPCs_Oligodendrocytes) <- "Condition"
Idents(Endothelial_cells) <- "Condition"
Idents(pilot.integrated) <- "Condition"

#Fresh cells
GO_Astrocytes_RGL_Fresh_cells<-DEenrichRPlot(Astrocytes_RGL, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
              max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_nIPCs_Neuroblasts_Fresh_cells<-DEenrichRPlot(nIPCs_Neuroblasts, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                 max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_Neurons_Fresh_cells<-DEenrichRPlot(Neurons, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                    max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_Microglia_Fresh_cells<-DEenrichRPlot(Microglia, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                          max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_OPCs_Oligodendrocytes_Fresh_cells<-DEenrichRPlot(OPCs_Oligodendrocytes, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                            max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_Endothelial_cells_Fresh_cells<-DEenrichRPlot(Endothelial_cells, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                        max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
#Fixed cells
GO_Astrocytes_RGL_Fixed_cells<-DEenrichRPlot(Astrocytes_RGL, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                             max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_nIPCs_Neuroblasts_Fixed_cells<-DEenrichRPlot(nIPCs_Neuroblasts, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_Neurons_Fixed_cells<-DEenrichRPlot(Neurons, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                      max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_Microglia_Fixed_cells<-DEenrichRPlot(Microglia, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                        max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_OPCs_Oligodendrocytes_Fixed_cells<-DEenrichRPlot(OPCs_Oligodendrocytes, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                    max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_Endothelial_cells_Fixed_cells<-DEenrichRPlot(Endothelial_cells, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)


GO_pilot.integrated_Fresh_cells<-DEenrichRPlot(pilot.integrated, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)
GO_pilot.integrated_Fixed_cells<-DEenrichRPlot(pilot.integrated, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                               max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=TRUE)


#Plotting the graphs
Graph_GO_Astrocytes_RGL_Fresh_cells<-DEenrichRPlot(Astrocytes_RGL, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                 max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_nIPCs_Neuroblasts_Fresh_cells<-DEenrichRPlot(nIPCs_Neuroblasts, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                    max.genes=515, enrich.database= "GO_Biological_Process_2018",return.gene.list=FALSE)
Graph_GO_Neurons_Fresh_cells<-DEenrichRPlot(Neurons, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                          max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_Microglia_Fresh_cells<-DEenrichRPlot(Microglia, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                            max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_OPCs_Oligodendrocytes_Fresh_cells<-DEenrichRPlot(OPCs_Oligodendrocytes, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                        max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_Endothelial_cells_Fresh_cells<-DEenrichRPlot(Endothelial_cells, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                    max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)

Graph_GO_Astrocytes_RGL_Fixed_cells<-DEenrichRPlot(Astrocytes_RGL, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                       max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_nIPCs_Neuroblasts_Fixed_cells<-DEenrichRPlot(nIPCs_Neuroblasts, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                          max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_Neurons_Fixed_cells<-DEenrichRPlot(Neurons, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_Microglia_Fixed_cells<-DEenrichRPlot(Microglia, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                  max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_OPCs_Oligodendrocytes_Fixed_cells<-DEenrichRPlot(OPCs_Oligodendrocytes, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                              max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_Endothelial_cells_Fixed_cells<-DEenrichRPlot(Endothelial_cells, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                          max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)


Graph_GO_pilot.integrated_Fresh_cells<-DEenrichRPlot(pilot.integrated, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                               max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)
Graph_GO_pilot.integrated_Fixed_cells<-DEenrichRPlot(pilot.integrated, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                               max.genes=515, enrich.database= "GO_Biological_Process_2021",return.gene.list=FALSE)


Graph_GO_pilot.integrated_Fresh_cells
Graph_GO_pilot.integrated_Fixed_cells

Graph_GO_Astrocytes_RGL_Fresh_cells
Graph_GO_nIPCs_Neuroblasts_Fresh_cells
Graph_GO_Neurons_Fresh_cells
Graph_GO_Microglia_Fresh_cells
Graph_GO_OPCs_Oligodendrocytes_Fresh_cells
Graph_GO_Endothelial_cells_Fresh_cells

Graph_GO_Astrocytes_RGL_Fixed_cells
Graph_GO_nIPCs_Neuroblasts_Fixed_cells
Graph_GO_Neurons_Fixed_cells
Graph_GO_Microglia_Fixed_cells
Graph_GO_OPCs_Oligodendrocytes_Fixed_cells
Graph_GO_Endothelial_cells_Fixed_cells

#Exporting GO analysis into an excel file (pendent de fer)


time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(GO_pilot.integrated_Fresh_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                        paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_pilot.integrated_Fresh_cells", col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(GO_pilot.integrated_Fixed_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_pilot.integrated_Fixed_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Astrocytes_RGL_Fresh_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Astrocytes_RGL_Fresh_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Astrocytes_RGL_Fixed_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Astrocytes_RGL_Fixed_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_nIPCs_Neuroblasts_Fresh_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_nIPCs_Neuroblasts_Fresh_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_nIPCs_Neuroblasts_Fixed_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_nIPCs_Neuroblasts_Fixed_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Neurons_Fresh_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Neurons_Fresh_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Neurons_Fixed_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Neurons_Fixed_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Microglia_Fresh_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Microglia_Fresh_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Microglia_Fixed_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Microglia_Fixed_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_OPCs_Oligodendrocytes_Fresh_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_OPCs_Oligodendrocytes_Fresh_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_OPCs_Oligodendrocytes_Fixed_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_OPCs_Oligodendrocytes_Fixed_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Endothelial_cells_Fresh_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Endothelial_cells_Fresh_cells", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Endothelial_cells_Fixed_cells, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Endothelial_cells_Fixed_cells", col.names = TRUE, row.names = TRUE, append = TRUE)









# GO_Cellular_Component_2021

Graph_GO_pilot.integrated_Fresh_cells_CC<-DEenrichRPlot(pilot.integrated, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                     max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_pilot.integrated_Fixed_cells_CC<-DEenrichRPlot(pilot.integrated, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                     max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_pilot.integrated_Fresh_cells_CC
Graph_GO_pilot.integrated_Fixed_cells_CC


#Fresh cells
GO_Astrocytes_RGL_Fresh_cells_CC<-DEenrichRPlot(Astrocytes_RGL, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                             max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_nIPCs_Neuroblasts_Fresh_cells_CC<-DEenrichRPlot(nIPCs_Neuroblasts, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_Neurons_Fresh_cells_CC<-DEenrichRPlot(Neurons, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                      max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_Microglia_Fresh_cells_CC<-DEenrichRPlot(Microglia, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                        max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_OPCs_Oligodendrocytes_Fresh_cells_CC<-DEenrichRPlot(OPCs_Oligodendrocytes, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                    max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_Endothelial_cells_Fresh_cells_CC<-DEenrichRPlot(Endothelial_cells, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
#Fixed cells
GO_Astrocytes_RGL_Fixed_cells_CC<-DEenrichRPlot(Astrocytes_RGL, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                             max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_nIPCs_Neuroblasts_Fixed_cells_CC<-DEenrichRPlot(nIPCs_Neuroblasts, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_Neurons_Fixed_cells_CC<-DEenrichRPlot(Neurons, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                      max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_Microglia_Fixed_cells_CC<-DEenrichRPlot(Microglia, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                        max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_OPCs_Oligodendrocytes_Fixed_cells_CC<-DEenrichRPlot(OPCs_Oligodendrocytes, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                    max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_Endothelial_cells_Fixed_cells_CC<-DEenrichRPlot(Endothelial_cells, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)


GO_pilot.integrated_Fresh_cells_CC<-DEenrichRPlot(pilot.integrated, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                               max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)
GO_pilot.integrated_Fixed_cells_CC<-DEenrichRPlot(pilot.integrated, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                               max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=TRUE)


#Plotting the graphs
Graph_GO_Astrocytes_RGL_Fresh_cells_CC<-DEenrichRPlot(Astrocytes_RGL, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                   max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_nIPCs_Neuroblasts_Fresh_cells_CC<-DEenrichRPlot(nIPCs_Neuroblasts, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                      max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_Neurons_Fresh_cells_CC<-DEenrichRPlot(Neurons, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                            max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_Microglia_Fresh_cells_CC<-DEenrichRPlot(Microglia, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                              max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_OPCs_Oligodendrocytes_Fresh_cells_CC<-DEenrichRPlot(OPCs_Oligodendrocytes, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                          max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_Endothelial_cells_Fresh_cells_CC<-DEenrichRPlot(Endothelial_cells, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                      max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)

Graph_GO_Astrocytes_RGL_Fixed_cells_CC<-DEenrichRPlot(Astrocytes_RGL, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                   max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_nIPCs_Neuroblasts_Fixed_cells_CC<-DEenrichRPlot(nIPCs_Neuroblasts, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                      max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_Neurons_Fixed_cells_CC<-DEenrichRPlot(Neurons, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                            max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_Microglia_Fixed_cells_CC<-DEenrichRPlot(Microglia, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                              max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_OPCs_Oligodendrocytes_Fixed_cells_CC<-DEenrichRPlot(OPCs_Oligodendrocytes, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                          max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)
Graph_GO_Endothelial_cells_Fixed_cells_CC<-DEenrichRPlot(Endothelial_cells, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                      max.genes=515, enrich.database= "GO_Cellular_Component_2021",return.gene.list=FALSE)

Graph_GO_pilot.integrated_Fresh_cells_CC
Graph_GO_pilot.integrated_Fixed_cells_CC

Graph_GO_Astrocytes_RGL_Fresh_cells_CC
Graph_GO_nIPCs_Neuroblasts_Fresh_cells_CC
Graph_GO_Neurons_Fresh_cells_CC
Graph_GO_Microglia_Fresh_cells_CC
Graph_GO_OPCs_Oligodendrocytes_Fresh_cells_CC
Graph_GO_Endothelial_cells_Fresh_cells_CC

Graph_GO_Astrocytes_RGL_Fixed_cells_CC
Graph_GO_nIPCs_Neuroblasts_Fixed_cells_CC
Graph_GO_Neurons_Fixed_cells_CC
Graph_GO_Microglia_Fixed_cells_CC
Graph_GO_OPCs_Oligodendrocytes_Fixed_cells_CC
Graph_GO_Endothelial_cells_Fixed_cells_CC

#Exporting GO analysis into an excel file 


time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(GO_pilot.integrated_Fresh_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_pilot.integrated_Fresh_cells_CC", col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(GO_pilot.integrated_Fixed_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                          paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_pilot.integrated_Fixed_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Astrocytes_RGL_Fresh_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                        paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Astrocytes_RGL_Fresh_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Astrocytes_RGL_Fixed_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                        paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Astrocytes_RGL_Fixed_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_nIPCs_Neuroblasts_Fresh_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                           paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_nIPCs_Neuroblasts_Fresh_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_nIPCs_Neuroblasts_Fixed_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                           paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_nIPCs_Neuroblasts_Fixed_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Neurons_Fresh_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                 paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Neurons_Fresh_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Neurons_Fixed_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                 paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Neurons_Fixed_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Microglia_Fresh_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                   paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Microglia_Fresh_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Microglia_Fixed_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                   paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Microglia_Fixed_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_OPCs_Oligodendrocytes_Fresh_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                               paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_OPCs_Oligodendrocytes_Fresh_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_OPCs_Oligodendrocytes_Fixed_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                               paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_OPCs_Oligodendrocytes_Fixed_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Endothelial_cells_Fresh_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                           paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Endothelial_cells_Fresh_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(GO_Endothelial_cells_Fixed_cells_CC, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/6 Gene expression comparations/", 
                                                           paste(time.name), "GO_pilot_Cellular_Component_2021 6 Gene expression comparations.xlsx"))
           , sheetName = "GO_Endothelial_cells_Fixed_cells_CC", col.names = TRUE, row.names = TRUE, append = TRUE)



# GO_Molecular_Function_2021

Graph_GO_pilot.integrated_Fresh_cells_MF<-DEenrichRPlot(pilot.integrated, ident.1 = "Fresh_cells", ident.2 = "Fixed_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                        max.genes=515, enrich.database= "GO_Molecular_Function_2021",return.gene.list=FALSE)
Graph_GO_pilot.integrated_Fixed_cells_MF<-DEenrichRPlot(pilot.integrated, ident.1 ="Fixed_cells", ident.2 = "Fresh_cells",test.use="wilcox", logfc.threshold=1, min.pct = 0.5, 
                                                        max.genes=515, enrich.database= "GO_Molecular_Function_2021",return.gene.list=FALSE)
Graph_GO_pilot.integrated_Fresh_cells_MF
Graph_GO_pilot.integrated_Fixed_cells_MF








#PCA analysis
#To see if there is variance between the two condition groups  PCA and check the first several PCs. (Wang 2021 )
#once integrated we need to re-assign the order of the conditions

pilot.integrated$Condition <- factor(pilot.integrated$Condition, levels=c("Fresh_cells","Fixed_cells"))
Idents(pilot.integrated) <- "Condition"

#Perform linear dimensional reduction
pilot.integrated <- RunPCA(pilot.integrated, features = VariableFeatures(object = pilot.integrated))
# Examine and visualize PCA results a few different ways
ElbowPlot(pilot.integrated)
P1<-DimPlot(pilot.integrated,dims = c(1, 2), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                       axis.text=element_text(size=6),
                                                                       axis.title.x = element_text(size =10),
                                                                       axis.title.y = element_text(size =10))
P2<-DimPlot(pilot.integrated,dims = c(1, 3), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                       axis.text=element_text(size=6),
                                                                       axis.title.x = element_text(size =10),
                                                                       axis.title.y = element_text(size =10))
P3<-DimPlot(pilot.integrated,dims = c(1, 4), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                       axis.text=element_text(size=6),
                                                                       axis.title.x = element_text(size =10),
                                                                       axis.title.y = element_text(size =10))
P4<-DimPlot(pilot.integrated,dims = c(1, 5), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                       axis.text=element_text(size=6),
                                                                       axis.title.x = element_text(size =10),
                                                                       axis.title.y = element_text(size =10))
P5<-DimPlot(pilot.integrated,dims = c(1, 6), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                       axis.text=element_text(size=6),
                                                                       axis.title.x = element_text(size =10),
                                                                       axis.title.y = element_text(size =10))
P6<-DimPlot(pilot.integrated,dims = c(1, 7), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                       axis.text=element_text(size=6),
                                                                       axis.title.x = element_text(size =10),
                                                                       axis.title.y = element_text(size =10))
P7<-DimPlot(pilot.integrated,dims = c(1, 8), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                       axis.text=element_text(size=6),
                                                                       axis.title.x = element_text(size =10),
                                                                       axis.title.y = element_text(size =10))
P8<-DimPlot(pilot.integrated,dims = c(1, 9), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                       axis.text=element_text(size=6),
                                                                       axis.title.x = element_text(size =10),
                                                                       axis.title.y = element_text(size =10))
P9<-DimPlot(pilot.integrated,dims = c(1, 10), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                        axis.text=element_text(size=6),
                                                                        axis.title.x = element_text(size =10),
                                                                        axis.title.y = element_text(size =10))
P10<-DimPlot(pilot.integrated,dims = c(1, 11), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                         axis.text=element_text(size=6),
                                                                         axis.title.x = element_text(size =10),
                                                                         axis.title.y = element_text(size =10))
P11<-DimPlot(pilot.integrated,dims = c(3, 2), reduction = "pca", pt.size=0.3)+ theme(legend.position = "none", 
                                                                        axis.text=element_text(size=6),
                                                                        axis.title.x = element_text(size =7),
                                                                        axis.title.y = element_text(size =7))
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10, ncol = 2, nrow = 5)

