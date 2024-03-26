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
#Initial name of the file "RDS file without ERCC + AG Urbach thresholds + separate pre-processing"
# Load the pilot dataset
pilot.data<-readRDS("F:/Marta/Getting started with scRNAseq analysis/Pilot results/scd-22aeac/diagnostics/2021-10-11-JEN-AU.rds",refhook = NULL)

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

#Different quality check graphs

#COMPARATION GRAPHS
VlnPlot(pilot, features = c("nFeature_RNA"), sort = "Condition", y.max=7000)
VlnPlot(pilot, features = c("nCount_RNA"), sort = "Condition", y.max=46000)
VlnPlot(pilot, features = c("percent_mt"), sort = "Condition", y.max=100)
VlnPlot(pilot, features = c("log10_nFeature_RNA"), log = TRUE, sort = "Condition", y.max = 5.2)
VlnPlot(pilot, features = c("log10_nCount_RNA"), log = TRUE, sort = "Condition", y.max = 6)

plot1 <- FeatureScatter(pilot, feature1 = "nCount_RNA", feature2 = "percent_mt", group.by="Condition")
plot2 <- FeatureScatter(pilot, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Condition")
plot1 + xlim (c(0,46000)) + ylim(c(0, 100))
plot2 + geom_abline(intercept = 0, slope = 0.25) + xlim (c(0,46000)) + ylim(c(0, 7000))
#GRAPHS AFTER FILTERING WIITH ADJUSTED AXIS

VlnPlot(pilot, features = c("nFeature_RNA", "nCount_RNA" ,"percent_mt"), ncol = 3, sort = "Condition")
VlnPlot(pilot, features = c("log10_nFeature_RNA", "log10_nCount_RNA"),ncol = 2, log = TRUE, sort = "Condition")

plot3 <- FeatureScatter(pilot, feature1 = "nCount_RNA", feature2 = "percent_mt", group.by="Condition")
plot4 <- FeatureScatter(pilot, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Condition")
plot3 
plot4 + geom_abline(intercept = 0, slope = 0.25)

#Histograms 
#FIXED AXIS
hist(pilot$log10_nCount_RNA, breaks=500, xlim = c(0,5))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fresh_cells$log10_nCount_RNA, breaks=500, xlim = c(0,5))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fixed_cells$log10_nCount_RNA,breaks=500, xlim = c(0,5))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
#ADAPTED AXIS
hist(pilot$log10_nCount_RNA, breaks=500)
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fresh_cells$log10_nCount_RNA,breaks=500)
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fixed_cells$log10_nCount_RNA,breaks=500)
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")

#FIXED AXIS
hist(pilot$percent_mt, breaks=500, xlim = c(0,100))
hist(Fresh_cells$percent_mt,breaks=500, xlim = c(0,100))
hist(Fixed_cells$percent_mt,breaks=500, xlim = c(0,100))
#ADAPTED AXIS
hist(pilot$percent_mt, breaks=500)
hist(Fresh_cells$percent_mt,breaks=500)
hist(Fixed_cells$percent_mt ,breaks=500)

#FIXED AXIS
hist(pilot$log10_nFeature_RNA, breaks=500, xlim = c(0,4))
abline(v=2.69897, col="blue")
hist(Fresh_cells$log10_nFeature_RNA,breaks=500, xlim = c(0,4))
abline(v=2.69897, col="blue")
hist(Fixed_cells$log10_nFeature_RNA ,breaks=500, xlim = c(0,4))
abline(v=2.69897, col="blue")
#ADAPTED AXIS
hist(pilot$log10_nFeature_RNA, breaks=500)
abline(v=2.69897, col="blue")
hist(Fresh_cells$log10_nFeature_RNA,breaks=500)
abline(v=2.69897, col="blue")
hist(Fixed_cells$log10_nFeature_RNA ,breaks=500)
abline(v=2.69897, col="blue")

#FIXED AXIS
hist(pilot$total_ercc_reads, breaks=500, xlim = c(0,2500))
abline(v=500, col="blue")
hist(Fresh_cells$total_ercc_reads,breaks=500, xlim = c(0,2500))
abline(v=500, col="blue")
hist(Fixed_cells$total_ercc_reads ,breaks=500, xlim = c(0,2500))
abline(v=500, col="blue")
#ADAPTED AXIS
hist(pilot$total_ercc_reads, breaks=500)
abline(v=500, col="blue")
hist(Fresh_cells$total_ercc_reads,breaks=500)
abline(v=500, col="blue")
hist(Fixed_cells$total_ercc_reads,breaks=500)
abline(v=500, col="blue")

#FIXED AXIS
hist(pilot$ratio_mol_gen, breaks=500, xlim = c(0,55))
abline(v=1.2, col="blue")
hist(Fresh_cells$ratio_mol_gen,breaks=500, xlim = c(0,55))
abline(v=1.2, col="blue")
hist(Fixed_cells$ratio_mol_gen,breaks=500, xlim = c(0,55))
abline(v=1.2, col="blue")
#ADAPTED AXIS
hist(pilot$ratio_mol_gen, breaks=500)
abline(v=1.2, col="blue")
hist(Fresh_cells$ratio_mol_gen,breaks=500)
abline(v=1.2, col="blue")
hist(Fixed_cells$ratio_mol_gen,breaks=500)
abline(v=1.2, col="blue")

plot5<-FeatureScatter(pilot, feature1 = "log10_nCount_RNA", feature2 = "log_Var_expression_per_cell", group.by="Condition", pt.size=2)
plot5 


#Normalisation
Fresh_cells<-SCTransform(Fresh_cells,new.assay.name="SCT", variable.features.n=2000)
Fresh_cells

Fixed_cells<-SCTransform(Fixed_cells,new.assay.name="SCT", variable.features.n=2000)
Fixed_cells

#For the scatter plot (var-mean)

Fresh_cells_nCount_asmatrix<-as.matrix(Fresh_cells$nCount_SCT)
Mean_exp_Fresh_SCT<-colMeans(Fresh_cells_nCount_asmatrix)
Fresh_cells[["log_Mean_exp_Fresh_SCT"]] <- log10(Mean_exp_Fresh_SCT)
Var_exp_Fresh_SCT<-colVars(Fresh_cells_nCount_asmatrix)
Fresh_cells[["log_Var_exp_Fresh_SCT"]] <- log10(Var_exp_Fresh_SCT)

Fixed_cells_nCount_asmatrix<-as.matrix(Fixed_cells$nCount_SCT)
Mean_exp_Fixed_SCT<-colMeans(Fixed_cells_nCount_asmatrix)
Fixed_cells[["log_Mean_exp_Fixed_SCT"]] <-log10(Mean_exp_Fixed_SCT)
Var_exp_Fixed_SCT<-colVars(Fixed_cells_nCount_asmatrix)
Fixed_cells[["log_Var_exp_Fixed_SCT"]] <- log10(Var_exp_Fixed_SCT)



#Scatter plots
plot6.1<-FeatureScatter(Fresh_cells, feature1 = "log10_nCount_RNA", feature2 = "log_Var_exp_Fresh_SCT", group.by="Condition", pt.size=2)
plot6.2<-FeatureScatter(Fixed_cells, feature1 = "log10_nCount_RNA", feature2 = "log_Var_exp_Fixed_SCT", group.by="Condition", pt.size=2, 
                        cols = "#00B8E7")
plot6.1
plot6.2


# Identify the 10 most highly variable genes
top10_Fresh_cells <- head(VariableFeatures(Fresh_cells, selection.method= NULL), 10)
top10_Fixed_cells <- head(VariableFeatures(Fixed_cells, selection.method= NULL), 10)
# plot variable features with and without labels
plot7.1 <- VariableFeaturePlot(Fresh_cells, selection.method= "sctransform")
plot7.2 <- LabelPoints(plot = plot7.1, points = top10_Fresh_cells, repel = TRUE) + ggtitle("Fresh_cells")
plot7.2

plot8.1 <- VariableFeaturePlot(Fixed_cells, selection.method= "sctransform")
plot8.2 <- LabelPoints(plot = plot8.1, points = top10_Fixed_cells, repel = TRUE)+ ggtitle("Fixed_cells")
plot8.2

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

#Perform linear dimensional reduction
pilot.integrated <- RunPCA(pilot.integrated, features = VariableFeatures(object = pilot.integrated))
# Examine and visualize PCA results a few different ways
DimPlot(pilot.integrated, reduction = "pca")
#Determine the ‘dimensionality’ of the dataset
#pilot.integrated <- JackStraw(pilot.integrated, num.replicate = 100)
#pilot.integrated <- ScoreJackStraw(pilot.integrated, dims = 1:20)
#JackStrawPlot(pilot.integrated, dims = 1:20)
ElbowPlot(pilot.integrated)




#Cluster the cells (12 dims)
pilot.integrated <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=10)
pilot.integrated<- FindClusters(pilot.integrated, resolution = 0.5, algorithm=1)


#Run non-linear dimensional reduction (UMAP/tSNE)

#tSNE

pilot.integrated<-RunTSNE(object = pilot.integrated)

DimPlot(pilot.integrated, reduction = "tsne")

#UMAP
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:12)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pilot.integrated, reduction = "umap")


#Lists of differential expressed genes (clusters)
pilot.markers <- FindAllMarkers(object=pilot.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox")
pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =2, order_by = avg_log2FC)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(pilot.markers, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot with Seurat + AG Urbach thresholds + separate pre-processing/", paste(time.name), "pilot.markers with AG Urbach thresholds + separate pre-processing.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)

#Differential expressed genes heatmap (clusters)
pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =10, order_by = avg_log2FC) -> top10_log2FC
top10_log2FC_HM<-DoHeatmap(pilot.integrated, features = top10_log2FC$gene)
top10_log2FC_HM + labs(title="top10_log2FC
                       ")

#Differential expressed genes heatmap (clusters)
pilot.markers %>%
  group_by(cluster) %>%
  slice_min(n =10, order_by = p_val) -> top10_pval
top10_pval_HM<-DoHeatmap(pilot.integrated, features = top10_pval$gene)
top10_pval_HM + labs(title="top10_pval
                       ")

#After cluster annotation with Linnarsson's database
#Assigning cell type identity to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("Astrocytes +  RGL", "Neuroblast 1 + Neuroblast 2 + nIPC ", "Granule-immature + Granule-mature + Neuroblast 2 + Mossy cells", 
                     "Granule-immature + Granule-mature + Mossy cells","Microglia","OL", "OPC + NFOL", "Endothelial")
pilot.integrated@active.ident <- plyr::mapvalues(x = pilot.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)

#UMAP
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:12)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pilot.integrated, reduction = "umap", label=TRUE)

#Differential expressed genes heatmap (clusters)
pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =10, order_by = avg_log2FC) -> top10_log2FC
top10_log2FC_HM<-DoHeatmap(pilot, features = top10_log2FC$gene)
top10_log2FC_HM + labs(title="top10_log2FC
                       ")

#Differential expressed genes heatmap (clusters)
pilot.markers %>%
  group_by(cluster) %>%
  slice_min(n =10, order_by = p_val) -> top10_pval
top10_pval_HM<-DoHeatmap(pilot, features = top10_pval$gene)
top10_pval_HM + labs(title="top10_pval
                       ")

# How do I create a UMAP plot where cells are colored by replicate?  First, store the current
# identities in a new column of meta.data called CellType
pilot.integrated$CellType <- Idents(pilot.integrated)
# Next, switch the identity class of all cells to reflect replicate ID
Idents(pilot.integrated) <- "Condition"
DimPlot(pilot.integrated, reduction = "umap", pt.size=1.7)
Idents(pilot.integrated) <- "CellType"
DimPlot(pilot.integrated, reduction = "umap", pt.size=1.7, label=TRUE)

FeaturePlot(pilot.integrated, features = c("nCount_RNA", "log10_nCount_RNA", "total_mito_reads", "percent.mt","nFeature_RNA", "log10_nFeature_RNA"))

#stacked barplots to see Contition contribution to cluster formation
pilot.integrated$Condition <- factor(pilot.integrated$Condition, levels=c("Fresh_cells","Fixed_cells"))
dittoBarPlot(pilot.integrated,"Condition", group.by = "CellType",retain.factor.levels=TRUE, color.panel=c("#F8766D","#00BFC4"))

#Graphs for progress report 21.02.2022


VlnPlot(pilot, features= "nFeature_RNA", sort = "Condition")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)
VlnPlot(pilot, features = c("nCount_RNA"), sort = "Condition")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)
VlnPlot(pilot, features = c("percent_mt"), sort = "Condition")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)

#UMAP_n15
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15,reduction.name="UMAP_n15",n.neighbors = 15)
Idents(pilot.integrated) <- "Condition"
DimPlot(pilot.integrated, reduction = "UMAP_n15", pt.size=2.3)
Idents(pilot.integrated) <- "CellType"
DimPlot(pilot.integrated, reduction = "UMAP_n15", pt.size=2.3) 

FeaturePlot(pilot.integrated, reduction= "UMAP_n15", features = c("Hopx--chr5", "Lpar1--chr4", "Prom1--chr5"), pt.size=2.3)

---------------------- 
#not adapted to integrated dataset

#Launch an interactive FeaturePlot
FeaturePlot(pilot.integrated, reduction= "UMAP_n15", features = c("nFeature_RNA"),interactive=TRUE)

#Now I will try to modify the UMAP appearance to make the cluster appear a little bit closer
#UMAP_n20
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15,reduction.name="UMAP_n20",n.neighbors = 20)
DimPlot(pilot.integrated, reduction = "UMAP_n20", pt.size=1.7, label=TRUE)
#UMAP_n15
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15,reduction.name="UMAP_n15",n.neighbors = 15)
DimPlot(pilot.integrated, reduction = "UMAP_n15", pt.size=1.7, label=TRUE)
#UMAP_n10
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15,reduction.name="UMAP_n10",n.neighbors = 10)
DimPlot(pilot.integrated, reduction = "UMAP_n10", pt.size=1.7, label=TRUE)
#UMAP_n11
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15,reduction.name="UMAP_n11",n.neighbors = 11)
DimPlot(pilot.integrated, reduction = "UMAP_n11", pt.size=1.7, label=TRUE)
#UMAP_n50
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15,reduction.name="UMAP_n50",n.neighbors = 50)
DimPlot(pilot.integrated, reduction = "UMAP_n50", pt.size=1.7, label=TRUE)
#UMAP_n5
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:13,reduction.name="UMAP_n5",n.neighbors = 5)
DimPlot(pilot.integrated, reduction = "UMAP_n5", pt.size=1.7, label=TRUE)
#UMAP_n199
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15,reduction.name="UMAP_n199",n.neighbors = 199)
DimPlot(pilot.integrated, reduction = "UMAP_n199", pt.size=1.7, label=TRUE)
#UMAP_n100
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15,reduction.name="UMAP_n100",n.neighbors = 100)
DimPlot(pilot.integrated, reduction = "UMAP_n100", pt.size=1.7, label=TRUE)
#UMAP_n101
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:13,reduction.name="UMAP_n101",n.neighbors = 101)
DimPlot(pilot.integrated, reduction = "UMAP_n101", pt.size=1.7, label=TRUE)
#UMAP_n150
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15,reduction.name="UMAP_n150",n.neighbors = 150)
DimPlot(pilot.integrated, reduction = "UMAP_n150", pt.size=1.7, label=TRUE)

