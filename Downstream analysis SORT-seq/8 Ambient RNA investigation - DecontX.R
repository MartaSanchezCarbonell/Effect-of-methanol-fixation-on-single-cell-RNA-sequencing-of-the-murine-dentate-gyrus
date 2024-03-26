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
BiocManager::install("celda")
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
install.packages("SoupX")
install.packages("vioplot")

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
library(SoupX)
library(celda)
library(vioplot)
#Load the seurat objects to extract the metadata for the SoupChannel object
pilot_tod<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/pilot_tod.rds")
Fresh_cells_tod<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/Fresh_cells_tod.rds")
Fixed_cells_tod<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/Fixed_cells_tod.rds")

pilot.integrated<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/pilot.integrated.rds")
Fresh_cells.final<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/Fresh_cells.final.rds")
Fixed_cells.final<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/Fixed_cells.final.rds")

#Create separate matrices for Ambient analysis per separate

pilot_mat<-as.matrix(pilot_tod@assays[["RNA"]]@counts)
Fresh_cells_mat<-as.matrix(Fresh_cells_tod@assays[["RNA"]]@counts)
Fixed_cells_mat<-as.matrix(Fixed_cells_tod@assays[["RNA"]]@counts)

pilot.integrated.mat<-as.matrix(pilot.integrated@assays[["RNA"]]@counts)
Fresh_cells.filt.mat<-as.matrix(Fresh_cells.final@assays[["RNA"]]@counts)
Fixed_cells.filt.mat<-as.matrix(Fixed_cells.final@assays[["RNA"]]@counts)


#create SCE objects
sce.pilot.raw<-SingleCellExperiment(list(counts = pilot_mat))
sce.Fresh_cells.raw<-SingleCellExperiment(list(counts = Fresh_cells_mat))
sce.Fixed_cells.raw<-SingleCellExperiment(list(counts = Fixed_cells_mat))

sce.pilot<-SingleCellExperiment(list(counts = pilot.integrated.mat))
sce.Fresh_cells<-SingleCellExperiment(list(counts = Fresh_cells.filt.mat))
sce.Fixed_cells<-SingleCellExperiment(list(counts = Fixed_cells.filt.mat))


#Run DecontX

sce.pilot<-decontX(sce.pilot,z=pilot.integrated$CellType,background =sce.pilot.raw, varGenes = 2000 )
sce.Fresh_cells<-decontX(sce.Fresh_cells,z=Fresh_cells.final$CellType,background =sce.Fresh_cells.raw, varGenes = 2000 )
sce.Fixed_cells<-decontX(sce.Fixed_cells,z=Fixed_cells.final$CellType,background =sce.Fixed_cells.raw, varGenes = 2000 )

#Violin plots of the contamination

vioplot(sce.Fresh_cells@metadata[["decontX"]][["contamination"]],border = "black", 
        col = 2,
        pchMed = 16,           
        plotCentre = "line",side = "left")
#stripchart(sce.Fresh_cells@metadata[["decontX"]][["contamination"]],side = "left", method = "jitter", col = "black",vertical = TRUE, pch = 19, add = TRUE)

vioplot(sce.Fixed_cells@metadata[["decontX"]][["contamination"]],border = "black", 
        col = 5,
        pchMed = 16,           
        plotCentre = "line",side = "right",add = TRUE)
#stripchart(sce.Fresh_cells@metadata[["decontX"]][["contamination"]], method = "jitter", col = "black",vertical = TRUE, pch = 19, add = TRUE)


#STATISTICS

#Means
Fresh_cells_contamination_mean<-mean(sce.Fresh_cells@metadata[["decontX"]][["contamination"]])
Fresh_cells_contamination_mean
Fixed_cells_contamination_mean<-mean(sce.Fixed_cells@metadata[["decontX"]][["contamination"]])
Fixed_cells_contamination_mean

#Media
Fresh_cells_contamination_median<-median(sce.Fresh_cells@metadata[["decontX"]][["contamination"]])
Fresh_cells_contamination_median
Fixed_cells_contamination_median<-median(sce.Fixed_cells@metadata[["decontX"]][["contamination"]])
Fixed_cells_contamination_median

#standard error of the mean (SEM)

std <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

Fresh_cells_contamination_std<-std(sce.Fresh_cells@metadata[["decontX"]][["contamination"]])
Fresh_cells_contamination_std
Fixed_cells_contamination_std<-std(sce.Fixed_cells@metadata[["decontX"]][["contamination"]])
Fixed_cells_contamination_std

#Summary
Condition<-c("Fresh cells","Fixed cells")
Contamination_mean<-c(Fresh_cells_contamination_mean,Fixed_cells_contamination_mean)
Contamination_median<-c(Fresh_cells_contamination_median,Fixed_cells_contamination_median)
Contamination_std<-c(Fresh_cells_contamination_std,Fixed_cells_contamination_std)
Summary<-data.frame(Condition,Contamination_mean,Contamination_median,Contamination_std)
Summary


shapiro.test(sce.Fresh_cells@metadata[["decontX"]][["contamination"]])
shapiro.test(sce.Fixed_cells@metadata[["decontX"]][["contamination"]])

wilcox.test(sce.Fresh_cells@metadata[["decontX"]][["contamination"]], sce.Fixed_cells@metadata[["decontX"]][["contamination"]], paired=FALSE)


#Checking how mitocondrial % is affected by Decontex

Fresh_cells.seurat<-as.Seurat(sce.Fresh_cells,counts = "counts", data = NULL)
Fixed_cells.seurat<-as.Seurat(sce.Fixed_cells,counts = "counts", data = NULL)

Fresh_cells.seurat[["percent_mt"]] <- PercentageFeatureSet(Fresh_cells.seurat, pattern = "^mt")
Fixed_cells.seurat[["percent_mt"]] <- PercentageFeatureSet(Fixed_cells.seurat, pattern = "^mt")

VlnPlot(Fresh_cells.seurat, "percent_mt", cols=c("#F8766D"))+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5) + scale_y_continuous(limits = c(0,100))
VlnPlot(Fixed_cells.seurat, "percent_mt", cols=c("#00BFC4"))+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5) + scale_y_continuous(limits = c(0,100))

#MEANS

Fresh_cells.mt.mean<-mean(Fresh_cells.seurat@meta.data[["percent_mt"]])
Fixed_cells.mt.mean<-mean(Fixed_cells.seurat@meta.data[["percent_mt"]])


Fresh_cells.mt.mean
Fixed_cells.mt.mean




#Downstream analysis-DecontX





#Loading into Seurat
Corrected_Fresh= CreateSeuratObject(decontXcounts(sce.Fresh_cells))
Corrected_Fixed= CreateSeuratObject(decontXcounts(sce.Fixed_cells))

#Checking how mitocondrial % is affected by soup removal

Corrected_Fresh[["percent_mt"]] <- PercentageFeatureSet(Corrected_Fresh, pattern = "^mt")
Corrected_Fixed[["percent_mt"]] <- PercentageFeatureSet(Corrected_Fixed, pattern = "^mt")

P1<-VlnPlot(Corrected_Fresh, "percent_mt", cols=c("#F8766D"))+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5) + theme(legend.position="none")+scale_y_continuous(limits = c(0,70))
P2<-VlnPlot(Corrected_Fixed, "percent_mt", cols=c("#00BFC4"))+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)+ theme(legend.position="none")+scale_y_continuous(limits = c(0,70))
ggarrange(P1,P2 ,ncol = 2, nrow = 1)
#means
Fresh_cells.mt.mean<-mean(Corrected_Fresh@meta.data[["percent_mt"]])
Fixed_cells.mt.mean<-mean(Corrected_Fixed@meta.data[["percent_mt"]])


Fresh_cells.mt.mean
Fixed_cells.mt.mean

#Normalisation
Corrected_Fresh<-SCTransform(Corrected_Fresh,new.assay.name="SCT", variable.features.n=2000)
Corrected_Fresh

Corrected_Fixed<-SCTransform(Corrected_Fixed,new.assay.name="SCT", variable.features.n=2000)
Corrected_Fixed

#Integration 
#Create a list of Seurat objects 
corrected_pilot.list1<-list(Corrected_Fresh,Corrected_Fixed)
# select features for downstream integration
corrected_pilot.features <- SelectIntegrationFeatures(object.list = corrected_pilot.list1, nfeatures = 2000)
#run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated
corrected_pilot.list1<-PrepSCTIntegration(object.list = corrected_pilot.list1, anchor.features = corrected_pilot.features, 
                                          verbose = FALSE)
#identify anchors
corrected_pilot.anchors<-FindIntegrationAnchors(object.list = corrected_pilot.list1, normalization.method = "SCT", 
                                                anchor.features = corrected_pilot.features, verbose = FALSE)
#integrate the datasets
corrected_pilot.integrated<-IntegrateData(anchorset = corrected_pilot.anchors, normalization.method = "SCT", 
                                          verbose = FALSE, k.weight =64)

#once integrated we need to re-assign the order of the conditions
corrected_pilot.integrated$Condition <- factor(pilot.integrated$Condition, levels=c("Corrected_Fresh","Corrected_Fixed"))

#Perform linear dimensional reduction
corrected_pilot.integrated <- RunPCA(corrected_pilot.integrated, features = VariableFeatures(object = corrected_pilot.integrated))
ElbowPlot(corrected_pilot.integrated)
#Cluster the cells (12 dims)
corrected_pilot.integrated <- FindNeighbors(corrected_pilot.integrated, dims = 1:12, k.param=10)
corrected_pilot.integrated<- FindClusters(corrected_pilot.integrated, resolution = 0.5, algorithm=1)

#UMAP
corrected_pilot.integrated <- RunUMAP(corrected_pilot.integrated, dims = 1:15, n.neighbors = 15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(corrected_pilot.integrated, reduction = "umap", label=TRUE, pt.size=2.3)



#Lists of differential expressed genes (clusters)
Corrected_pilot.markers <- FindAllMarkers(object=corrected_pilot.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox", assay="RNA")
Corrected_pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =2, order_by = avg_log2FC)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Corrected_pilot.markers, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/", paste(time.name), "pilot.markers corrected with soup X.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)

#Assigning cell type identity to clusters
corrected_current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
corrected_new.cluster.ids <- c("Astrocytes +  RGL", "nIPCs + Neuroblasts", "Immature neurons", 
                               "Mature neurons","Microglia","Oligodendrocytes", "OPCs", "Endothelial")
corrected_pilot.integrated@active.ident <- plyr::mapvalues(x = corrected_pilot.integrated@active.ident, from = corrected_current.cluster.ids, to = corrected_new.cluster.ids)
corrected_pilot.integrated$CellType <- Idents(corrected_pilot.integrated)

#UMAP
corrected_pilot.integrated <- RunUMAP(corrected_pilot.integrated, dims = 1:15, n.neighbors = 15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(corrected_pilot.integrated, reduction = "umap", label=TRUE, pt.size=2.3)


#Differential expressed genes heatmap (clusters)
Corrected_pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =10, order_by = avg_log2FC) -> top10_log2FC
top10_log2FC_HM<-DoHeatmap(corrected_pilot.integrated, features = top10_log2FC$gene, label=FALSE)
top10_log2FC_HM + labs(title="top10_log2FC
                       ")

#Differential expressed genes heatmap (clusters)
Corrected_pilot.markers %>%
  group_by(cluster) %>%
  slice_min(n =10, order_by = p_val) -> top10_pval
top10_pval_HM<-DoHeatmap(corrected_pilot.integrated, features = top10_pval$gene, label=FALSE)
top10_pval_HM + labs(title="top10_pval
                       ")

#Differential expressed genes heatmap (clusters)
Corrected_pilot.markers %>%
  group_by(cluster) %>%
  slice_min(n =10, order_by = p_val_adj) -> top10_pval_adj
top10_pval_adj_HM<-DoHeatmap(corrected_pilot.integrated, features = top10_pval_adj$gene, label=FALSE)
top10_pval_adj_HM + labs(title="top10_pval_adj")
























































#Creation of the Seurat Objects to load the cluster information for DecontX


# Load the pilot dataset
pilot.data<-readRDS("F:/Marta/Getting started with scRNAseq analysis/Pilot results/scd-22aeac/diagnostics/2021-10-11-JEN-AU.rds",refhook = NULL)

#remove chromosome lables
row.names(pilot.data)<-gsub("__.....","",row.names(pilot.data))
row.names(pilot.data)<-gsub("__....","",row.names(pilot.data))
#Load the seurat objects to extract the metadata for the SoupChannel object

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
pilot_tod <- CreateSeuratObject(counts = pilot.data[-ercc,], project = "Pilot analysis", min.cells = 0, min.features = 0)
pilot_tod


pilot_tod$orig.ident <- factor(pilot_tod$orig.ident, levels=c("JEN-AU-s004","JEN-AU-s002"))

current.ids <- c("JEN-AU-s004","JEN-AU-s002")
new.ids <- c("Fresh_cells","Fixed_cells")
pilot_tod@active.ident <- plyr::mapvalues(x = pilot_tod@active.ident, from = current.ids, to = new.ids)
pilot_tod$Condition <- Idents(pilot_tod)

pilot_tod$Condition <- factor(pilot_tod$Condition, levels=c("Fresh_cells","Fixed_cells"))


#Pre-processing workflow

#Preparing all the extra columns in our surat Object
#QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pilot_tod[["Var_expression_per_cell"]] <- Var_expression_per_cell
pilot_tod[["Mean_expression_per_cell"]] <- Mean_expression_per_cell
pilot_tod[["log_Var_expression_per_cell"]] <- log10(pilot_tod$Var_expression_per_cell)
pilot_tod[["log_Mean_expression_per_cell"]] <- log10(pilot_tod$Mean_expression_per_cell)
pilot_tod[["total_mito_reads"]] <- total_mito_reads
pilot_tod[["total_ercc_reads"]] <- total_ercc_reads
pilot_tod[["percent_ERCC"]] <- percent.ERCC
pilot_tod[["percent_mt"]] <- PercentageFeatureSet(pilot_tod, pattern = "^mt")
pilot_tod[["percent_ERCC"]] <- PercentageFeatureSet(pilot_tod, pattern = "^ERCC")
pilot_tod[["log10_nCount_RNA"]] <- log10(pilot_tod$nCount_RNA)
pilot_tod[["log10_nFeature_RNA"]] <- log10(pilot_tod$nFeature_RNA)
pilot_tod[["percent_ERCC_after_removal"]] <- PercentageFeatureSet(pilot_tod, pattern = "^ERCC")
pilot_tod[["ratio_mol_gen"]]<-pilot_tod$nCount_RNA/pilot_tod$nFeature_RNA
Control_of_ERCC_removal<- summary(pilot_tod$percent_ERCC_after_removal)
Control_of_ERCC_removal


pilot_tod.list<-SplitObject(pilot_tod, split.by = "Condition")

Fresh_cells_tod<-pilot_tod.list$Fresh_cells
Fixed_cells_tod<-pilot_tod.list$Fixed_cells


#Saving tod objects

saveRDS(pilot_tod,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/pilot_tod.rds")
saveRDS(Fresh_cells_tod,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/Fresh_cells_tod.rds")
saveRDS(Fixed_cells_tod,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/Fixed_cells_tod.rds")





#filtering
Fresh_cells<-subset(Fresh_cells_tod, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
Fixed_cells<-subset(Fixed_cells_tod, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
pilot<-subset(pilot_tod, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
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
pilot.integrated$Condition <- factor(pilot.integrated$Condition, levels=c("Fresh_cells","Fixed_cells"))

#Perform linear dimensional reduction
pilot.integrated <- RunPCA(pilot.integrated, features = VariableFeatures(object = pilot.integrated))
#Cluster the cells (12 dims)
pilot.integrated <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=10)
pilot.integrated<- FindClusters(pilot.integrated, resolution = 0.5, algorithm=1)
#After cluster annotation with Linnarsson's database
#Assigning cell type identity to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("Astrocytes +  RGL", "nIPCs + Neuroblasts", "Immature neurons ", 
                     "Mature neurons","Microglia","Oligodendrocytes", "OPCs", "Endothelial")
pilot.integrated@active.ident <- plyr::mapvalues(x = pilot.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
pilot.integrated$CellType <- Idents(pilot.integrated)

#UMAP
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15, n.neighbors = 15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pilot.integrated, reduction = "umap", label=TRUE, pt.size=2.3)

pilot.final.list<-SplitObject(pilot.integrated, split.by = "Condition")

Fresh_cells.final<-pilot.final.list$Fresh_cells
Fixed_cells.final<-pilot.final.list$Fixed_cells

DimPlot(Fresh_cells.final, reduction = "umap", label=TRUE, pt.size=2.3)
DimPlot(Fixed_cells.final, reduction = "umap", label=TRUE, pt.size=2.3)


saveRDS(pilot.integrated,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/pilot.integrated.rds")
saveRDS(Fresh_cells.final,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/Fresh_cells.final.rds")
saveRDS(Fixed_cells.final,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/soupX objects/Fixed_cells.final.rds")

