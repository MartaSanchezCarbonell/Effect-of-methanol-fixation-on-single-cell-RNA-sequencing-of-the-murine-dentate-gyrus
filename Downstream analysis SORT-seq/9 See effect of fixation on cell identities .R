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
Idents(pilot) <- "Condition"

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
#VlnPlot(pilot, features = c("nFeature_RNA"), sort = "Condition", y.max=7000)
#VlnPlot(pilot, features = c("nCount_RNA"), sort = "Condition", y.max=46000)
#VlnPlot(pilot, features = c("percent_mt"), sort = "Condition", y.max=100)
#VlnPlot(pilot, features = c("log10_nFeature_RNA"), log = TRUE, sort = "Condition", y.max = 5.2)
#VlnPlot(pilot, features = c("log10_nCount_RNA"), log = TRUE, sort = "Condition", y.max = 6)

#plot1 <- FeatureScatter(pilot, feature1 = "nCount_RNA", feature2 = "percent_mt", group.by="Condition")
#plot2 <- FeatureScatter(pilot, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Condition")
#plot1 + xlim (c(0,46000)) + ylim(c(0, 100))
#plot2 + geom_abline(intercept = 0, slope = 0.25) + xlim (c(0,46000)) + ylim(c(0, 7000))
#GRAPHS AFTER FILTERING WIITH ADJUSTED AXIS

#VlnPlot(pilot, features = c("nFeature_RNA", "nCount_RNA" ,"percent_mt"), ncol = 3, sort = "Condition")
#VlnPlot(pilot, features = c("log10_nFeature_RNA", "log10_nCount_RNA"),ncol = 2, log = TRUE, sort = "Condition")

#plot3 <- FeatureScatter(pilot, feature1 = "nCount_RNA", feature2 = "percent_mt", group.by="Condition")
#plot4 <- FeatureScatter(pilot, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Condition")
#plot3 
#plot4 + geom_abline(intercept = 0, slope = 0.25)

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

#plot5<-FeatureScatter(pilot, feature1 = "log10_nCount_RNA", feature2 = "log_Var_expression_per_cell", group.by="Condition", pt.size=2)
#plot5 


#Normalisation
Fresh_cells<-SCTransform(Fresh_cells,new.assay.name="SCT", variable.features.n=2000)
Fresh_cells

Fixed_cells<-SCTransform(Fixed_cells,new.assay.name="SCT", variable.features.n=2000)
Fixed_cells

#For the scatter plot (var-mean)

#Fresh_cells_nCount_asmatrix<-as.matrix(Fresh_cells$nCount_SCT)
#Mean_exp_Fresh_SCT<-colMeans(Fresh_cells_nCount_asmatrix)
#Fresh_cells[["log_Mean_exp_Fresh_SCT"]] <- log10(Mean_exp_Fresh_SCT)
#Var_exp_Fresh_SCT<-colVars(Fresh_cells_nCount_asmatrix)
#Fresh_cells[["log_Var_exp_Fresh_SCT"]] <- log10(Var_exp_Fresh_SCT)

#Fixed_cells_nCount_asmatrix<-as.matrix(Fixed_cells$nCount_SCT)
#Mean_exp_Fixed_SCT<-colMeans(Fixed_cells_nCount_asmatrix)
#Fixed_cells[["log_Mean_exp_Fixed_SCT"]] <-log10(Mean_exp_Fixed_SCT)
#Var_exp_Fixed_SCT<-colVars(Fixed_cells_nCount_asmatrix)
#Fixed_cells[["log_Var_exp_Fixed_SCT"]] <- log10(Var_exp_Fixed_SCT)



#Scatter plots
#plot6.1<-FeatureScatter(Fresh_cells, feature1 = "log10_nCount_RNA", feature2 = "log_Var_exp_Fresh_SCT", group.by="Condition", pt.size=2)
#plot6.2<-FeatureScatter(Fixed_cells, feature1 = "log10_nCount_RNA", feature2 = "log_Var_exp_Fixed_SCT", group.by="Condition", pt.size=2, 
                        #cols = "#00B8E7")
#plot6.1
#plot6.2


# Identify the 10 most highly variable genes
#top10_Fresh_cells <- head(VariableFeatures(Fresh_cells, selection.method= NULL), 10)
#top10_Fixed_cells <- head(VariableFeatures(Fixed_cells, selection.method= NULL), 10)
# plot variable features with and without labels
#plot7.1 <- VariableFeaturePlot(Fresh_cells, selection.method= "sctransform")
#plot7.2 <- LabelPoints(plot = plot7.1, points = top10_Fresh_cells, repel = TRUE) + ggtitle("Fresh_cells")
#plot7.2

#plot8.1 <- VariableFeaturePlot(Fixed_cells, selection.method= "sctransform")
#plot8.2 <- LabelPoints(plot = plot8.1, points = top10_Fixed_cells, repel = TRUE)+ ggtitle("Fixed_cells")
#plot8.2

#Integration 
#Create a list of Seurat objects 
pilot.list<-list(Fresh_cells,Fixed_cells)
# select features for downstream integration
pilot.features <- SelectIntegrationFeatures(object.list = pilot.list, nfeatures = 2000)
#run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated
pilot.list<-PrepSCTIntegration(object.list = pilot.list, anchor.features = pilot.features, 
                                verbose = FALSE)
#identify anchors
pilot.anchors<-FindIntegrationAnchors(object.list = pilot.list, normalization.method = "SCT", 
                                      anchor.features = pilot.features, verbose = FALSE)
#integrate the datasets
pilot.integrated<-IntegrateData(anchorset = pilot.anchors, normalization.method = "SCT", 
                                verbose = FALSE, k.weight =64)

#once integrated we need to re-assign the order of the conditions
pilot.integrated$Condition <- factor(pilot.integrated$Condition, levels=c("Fresh_cells","Fixed_cells"))
Idents(pilot.integrated) <- "Condition"


#Perform linear dimensional reduction
pilot.integrated <- RunPCA(pilot.integrated, features = VariableFeatures(object = pilot.integrated))
# Examine and visualize PCA results a few different ways
DimPlot(pilot.integrated, reduction = "pca")
#Determine the ‘dimensionality’ of the dataset
#pilot.integrated <- JackStraw(pilot.integrated, num.replicate = 100)
#pilot.integrated <- ScoreJackStraw(pilot.integrated, dims = 1:20)
#JackStrawPlot(pilot.integrated, dims = 1:20)
ElbowPlot(pilot.integrated)









#To see effect of fixation in distinct subpopulations 
#Cluster analysis





#Cluster the cells (12 dims)
pilot.integrated <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=10)
pilot.integrated<- FindClusters(pilot.integrated, resolution = 0.5, algorithm=1)






#Cluster annotation











#Lists of differential expressed genes (clusters)
pilot.markers <- FindAllMarkers(object=pilot.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox",assay="RNA")
pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =2, order_by = avg_log2FC)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(pilot.markers, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/9 See effect of fixation on cell identities/", paste(time.name), "pilot.markers 9 See effect of fixation on cell identities.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)



#After cluster annotation with Linnarsson's database
#Assigning cell type identity to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("Astrocytes +  RGL", "nIPCs + Neuroblasts", "Immature neurons", 
                     "Mature neurons","Microglia","Oligodendrocytes", "OPCs", "Endothelial")
pilot.integrated@active.ident <- plyr::mapvalues(x = pilot.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
pilot.integrated$CellType <- Idents(pilot.integrated)

#UMAP
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15, n.neighbors = 15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pilot.integrated, reduction = "umap", label=FALSE, pt.size=2.3)

#Differential expressed genes heatmap (clusters)
pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =10, order_by = avg_log2FC) -> top10_log2FC
top10_log2FC_HM<-DoHeatmap(pilot.integrated, features = top10_log2FC$gene, label=FALSE)
top10_log2FC_HM + labs(title="top10_log2FC
                       ")

#Differential expressed genes heatmap (clusters)
pilot.markers %>%
  group_by(cluster) %>%
  slice_min(n =10, order_by = p_val) -> top10_pval
top10_pval_HM<-DoHeatmap(pilot.integrated, features = top10_pval$gene, label=FALSE, assay = )
top10_pval_HM + labs(title="top10_pval
                       ")





#Plot marker gene expression on UMAP (keep 2000 variable genes and change the active assay 
#to RNA when searching for genes) (Chen 2018) (Xperience 2022) Expression coloured based on normalized 
#expression levels. (Alles 2017) (Wang 2021 )


#I will check for the markers in our UMAP with all cells 
#and I will subset and cluster the astrocyte and RGL cluster to see if I can separate the RGLs   
 
DefaultAssay(pilot.integrated) <- "RNA"
FeaturePlot(pilot.integrated, features = c("Lpar1", "Nes", "Prom1", "Hes5", "Ascl1", "Tfap2c", "Vnn1", "Wnt8b", "2700094K13Rik", 
                                           "Hopx", "Emx1", "E330013P04Rik", "Btg2", "Frzb","Rhcg", 
                                           "Sox4", "Riiad1", "Fabp7"), pt.size=1)










pilot.integrated$CellType <- Idents(pilot.integrated)
Astrocytes_RGL<-subset(pilot.integrated, subset = CellType=="Astrocytes +  RGL")
FeaturePlot(Astrocytes_RGL, features = c("Lpar1", "Nes", "Prom1", "Hes5", "Ascl1", "Tfap2c", "Vnn1", "Wnt8b", "2700094K13Rik", 
                                           "Hopx", "Emx1", "E330013P04Rik", "Btg2", "Frzb", "Rhcg",
                                         "Sox4", "Riiad1", "Fabp7"), pt.size=2)

#Cluster analysisi from Astrocytes_RGL!!!
DefaultAssay(Astrocytes_RGL) <- "integrated"
Idents(Astrocytes_RGL) <- "Condition"
#Perform linear dimensional reduction
#Astrocytes_RGL <- RunPCA(Astrocytes_RGL, features = VariableFeatures(object = Astrocytes_RGL))
# Examine and visualize PCA results a few different ways
DimPlot(Astrocytes_RGL, reduction = "pca")
#Determine the ‘dimensionality’ of the dataset
ElbowPlot(Astrocytes_RGL)

#Cluster analysis

#Cluster the cells (12 dims)
Astrocytes_RGL <- FindNeighbors(Astrocytes_RGL, dims = 1:12, k.param=10, graph.name = "k10")
Astrocytes_RGL<- FindClusters(Astrocytes_RGL, resolution = 0.5, algorithm=1, graph.name = "k10")

#UMAP
Astrocytes_RGL <- RunUMAP(Astrocytes_RGL, dims = 1:12, n.neighbors = 12)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(Astrocytes_RGL, reduction = "umap", label=TRUE, pt.size=2.3)


#Lists of differential expressed genes (clusters)
Astrocytes_RGL.markers <- FindAllMarkers(object=Astrocytes_RGL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox")
Astrocytes_RGL.markers %>%
  group_by(cluster) %>%
  slice_max(n =2, order_by = avg_log2FC)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Astrocytes_RGL.markers, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/9 See effect of fixation on cell identities/", paste(time.name), "Astrocytes_RGL.markers 9 See effect of fixation on cell identities.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)


#decreasing the K

#Cluster the cells (12 dims)
Astrocytes_RGL <- FindNeighbors(Astrocytes_RGL, dims = 1:12, k.param=4, graph.name = "k4")
Astrocytes_RGL<- FindClusters(Astrocytes_RGL, resolution = 0.5, algorithm=1, graph.name = "k4")


#UMAP
Astrocytes_RGL <- RunUMAP(Astrocytes_RGL, dims = 1:12, n.neighbors = 3)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(Astrocytes_RGL, reduction = "umap", label=TRUE, pt.size=2.3)

#Lists of differential expressed genes (clusters)
Astrocytes_RGL.markers_k4 <- FindAllMarkers(object=Astrocytes_RGL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox", assay="RNA")
Astrocytes_RGL.markers_k4 %>%
  group_by(cluster) %>%
  slice_max(n =2, order_by = avg_log2FC)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Astrocytes_RGL.markers_k4, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/9 See effect of fixation on cell identities/", paste(time.name), "Astrocytes_RGL.markers_k4 RNA assay 9 See effect of fixation on cell identities.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)









DefaultAssay(Astrocytes_RGL) <- "RNA"
FeaturePlot(Astrocytes_RGL, features = c("Lpar1", "Nes", "Prom1", "Hes5", "Ascl1", "Tfap2c", "Vnn1", "Wnt8b", "2700094K13Rik", 
                                         "Hopx", "Emx1", "E330013P04Rik", "Btg2", "Frzb", "Rhcg",
                                         "Sox4", "Riiad1", "Fabp7"), pt.size=2)


A1<-DimPlot(Astrocytes_RGL, reduction = "umap", label=TRUE, pt.size=2.3)

A2<-DimPlot(Astrocytes_RGL, reduction = "umap", label=FALSE, pt.size=2.3, group.by = "Condition")

A1+A2

#Row counts
DoHeatmap(Astrocytes_RGL, features = c("Lpar1", "Nes", "Prom1", "Hes5", "Ascl1", "Tfap2c", "Vnn1", "Wnt8b", "2700094K13Rik", 
                                       "Hopx", "Emx1", "E330013P04Rik", "Btg2", "Frzb", "Rhcg",
                                       "Sox4", "Riiad1", "Fabp7"), label=TRUE, assay = "RNA", slot = "counts")
#Corrected counts
DoHeatmap(Astrocytes_RGL, features = c("Lpar1", "Nes", "Prom1", "Hes5", "Ascl1", "Tfap2c", "Vnn1", "Wnt8b", "2700094K13Rik", 
                                       "Hopx", "Emx1", "E330013P04Rik", "Btg2", "Frzb", "Rhcg",
                                       "Sox4", "Riiad1", "Fabp7"), label=TRUE, assay = "SCT", slot = "counts")
# log2 counts
DoHeatmap(Astrocytes_RGL, features = c("Lpar1", "Nes", "Prom1", "Hes5", "Ascl1", "Tfap2c", "Vnn1", "Wnt8b", "2700094K13Rik", 
                                       "Hopx", "Emx1", "E330013P04Rik", "Btg2", "Frzb", "Rhcg",
                                       "Sox4", "Riiad1", "Fabp7"), label=TRUE, assay = "SCT", slot = "data")
#Pearson residuals
DoHeatmap(Astrocytes_RGL, features = c("Lpar1", "Nes", "Prom1", "Hes5", "Ascl1", "Tfap2c", "Vnn1", "Wnt8b", "2700094K13Rik", 
                                       "Hopx", "Emx1", "E330013P04Rik", "Btg2", "Frzb", "Rhcg",
                                       "Sox4", "Riiad1", "Fabp7"), label=TRUE, assay = "integrated", slot = "data")









#To see if methanol fixation maintain cell composition  Cell type composition graphs 
#(pie graphs with % of the different cell identities or graph of cluster contribution that we already have) (Denisenko 2020)













#UMAP plot where cells are colored by replicate?  
# First, store the current identities in a new column of meta.data called CellType
pilot.integrated$CellType <- Idents(pilot.integrated)
# Next, switch the identity class of all cells to reflect replicate ID
Idents(pilot.integrated) <- "Condition"
DimPlot(pilot.integrated, reduction = "umap", pt.size=2.3)
Idents(pilot.integrated) <- "CellType"
DimPlot(pilot.integrated, reduction = "umap", pt.size=2.3, split.by = "Condition")+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))




#stacked barplots to see Condition contribution to cluster formation
Idents(pilot.integrated) <- "Condition"
dittoBarPlot(pilot.integrated,var="Condition", group.by = "CellType",retain.factor.levels=TRUE, color.panel=c("#F8766D","#00BFC4"))

#pie graphs with % of the different cell identities
## extract meta data
md <- pilot.integrated@meta.data %>% as.data.table
# the resulting md object has one "row" per cell

## count the number of cells per unique combinations of "Condition" and "CellType" ad plot the pie chart
Proportion_CellType_Condition<-md[, .N, by = c("Condition", "CellType")]
Proportion_Fresh<-subset(Proportion_CellType_Condition, Condition=="Fresh_cells")
bp_Fresh<- ggplot(Proportion_Fresh, aes(x="", y=N, fill=CellType))+
  geom_bar(width = 1, stat = "identity")
bp_Fresh
pie_Fresh <- bp_Fresh + coord_polar("y", start=0) + ggtitle("Proportion_Fresh")+ geom_text(aes(label=CellType),position = position_stack(vjust = 0.5), size=3.5)
pie_Fresh

Proportion_Fixed<-subset(Proportion_CellType_Condition, Condition=="Fixed_cells")
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
write.xlsx(Proportion_CellType_Condition_Percentage, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/9 See effect of fixation on cell identities/", paste(time.name), "Proportion_CellType_Condition with Percentage.xlsx"))
           , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)


#Plots of with the typical cell type markers ONLY WITH FIXED CELLS!!!


#UMAP
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15, n.neighbors = 15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pilot.integrated, reduction = "umap", label=FALSE, pt.size=2.3)

pilot.list.int<-SplitObject(pilot.integrated, split.by = "Condition")

Fixed_cells.i<-pilot.list.int$Fixed_cells



DefaultAssay(Fixed_cells.i)<-"SCT"
#Granule mature 
FeaturePlot(Fixed_cells.i, features= c("Prox1","Plk5", "Ntng1", "Rbfox3", "Calb1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Prox1","Plk5", "Ntng1", "Rbfox3", "Calb1"), slot ="counts")
#Granule mature (round 2)
VlnPlot(Fixed_cells.i,ncol=6, features= c("Rgs4","Trpc6","Bc030499","Ntng1","Rasgrp1","Nos1","Snhg11","Smad3","Adarb2","Ipcef1","Kcnip3","Tekt5","Gm12216","Kart2","Hlf","Plekha2","Icam4","Mfsd4","Tanc1","Cstad"), slot ="counts")

#Granule immature
FeaturePlot(Fixed_cells.i, features= c("Ncam1","Prox1","Dcx", "Fxyd7","Rbfox3", "Calb2" ), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Ncam1","Prox1","Dcx", "Fxyd7","Rbfox3", "Calb2" ), slot ="counts")
#Granule immature (round 2)
VlnPlot(Fixed_cells.i,ncol=6, features= c("Camk4", "Stc1", "Sdk2", "Rasd1", "Fxyd7", "Il16", "Dsp", "Rspo3", "Prickle2", "Smim3", "Mcm6", "Rasgrf2", "Tenm1", "Scl4a4", "Gda", "Ptprj", "Rbm24","Omg", "Dnajc27"), slot ="counts")

#Neuroblast 2
FeaturePlot(Fixed_cells.i, features= c("Prox1","Dcx","Gal", "Sox11","Neurod1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Prox1","Dcx","Gal", "Sox11","Neurod1"), slot ="counts")
#Neuroblast1
FeaturePlot(Fixed_cells.i, features= c("Prox1","Dcx","Eomes", "Igfbpl1", "Calb2", "Neurog2", "Sox11", "Neurod4", "Tac2", "Neurod1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Prox1","Dcx","Eomes", "Igfbpl1", "Calb2", "Neurog2", "Sox11", "Neurod4", "Tac2", "Neurod1"), slot ="counts")
#nIPC 
FeaturePlot(Fixed_cells.i, features= c("Prox1","Dcx","Ascl", "Cdk1", "Tfap2c", "Eomes", "Igfbpl1","Neurog2", "Top2a", "Neurod4", "Aurkb", "Neurod1", "Ccnd2"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Prox1","Dcx","Ascl", "Cdk1", "Tfap2c", "Eomes", "Igfbpl1","Neurog2", "Top2a", "Neurod4", "Aurkb", "Neurod1", "Ccnd2"), slot ="counts")
#RGL
FeaturePlot(Fixed_cells.i, features= c("Tfap2c", "Sox9", "Lpar1", "Ascl1","Sox2", "Nes", "Vnn1","Hes5", "Rhcg", "Scl1a3", "Wnt8b", "Aldoc", "Apoe", "Id4", "Hopx", "GFAP"), slot ="counts")
VlnPlot(Fixed_cells.i, features= c("Tfap2c", "Sox9", "Lpar1", "Ascl1","Sox2", "Nes", "Vnn1","Hes5", "Rhcg", "Scl1a3", "Wnt8b", "Aldoc", "Apoe", "Id4", "Hopx", "GFAP"), slot ="counts")
#Astrocytes
FeaturePlot(Fixed_cells.i, features= c("S100b", "Fzd2", "Sox9","Sox2","Hes5","Scl1a3"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("S100b", "Fzd2", "Sox9","Sox2","Hes5","Scl1a3"), slot ="counts")
#OPC 
FeaturePlot(Fixed_cells.i, features= c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"), slot ="counts")
#Oligo 
FeaturePlot(Fixed_cells.i, features= c("Lpar1","Mbp", "Plp1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Lpar1","Mbp", "Plp1"), slot ="counts")
#Endothelial cells
FeaturePlot(Fixed_cells.i, features= c("Nes", "Prom1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, ncol=6, features= c("Nes", "Prom1", "Efna1", "Tm4sf1", "Pcp4l1","Cd93","Igfbp7","Hspb1","Flt1","Esam","Eltd1","Abcb1a","Cldn5"), slot ="counts")
#Microglia
FeaturePlot(Fixed_cells.i, features= c("Csf1r" , "Cx3cr1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Csf1r" , "Cx3cr1"), slot ="counts")


#Selected markers Heatmap

DoHeatmap(Fixed_cells.i, features = c("Rasgrp","Smad3","Kcnip3","Mfsd4","Tanc1","Dsp","Rasgrf2","Gda","Rbm24"), slot ="scale.data")


#To be able to see the markers and not have a to wide range of expression we need to plot the log2(1+counts) expression
#(if we had 1 cell with an expression of 40 we were no able to see cells with expression of 15)

#Create a new assay to store log2(1+counts) information
#Save the linear normalized (LN) counts matrix
LN_counts_Fixed_cells.i<-as.data.frame(Fixed_cells.i@assays[["SCT"]]@counts)
#add +1 to each element of the matrix
LN_counts_Fixed_cells.i_1<-LN_counts_Fixed_cells.i+1
#check that it worked
LN_counts_Fixed_cells.i[1:10,1:4]
LN_counts_Fixed_cells.i_1[1:10,1:4]
#do the log2 from the counts+1 matrix
LN_counts_Fixed_cells.i_1_log2<-log2(LN_counts_Fixed_cells.i_1)
LN_counts_Fixed_cells.i_1_log2[1:10,1:4]
#Create a new assay to store log2(1+counts) information
log2_assay <- CreateAssayObject(counts = LN_counts_Fixed_cells.i_1_log2)

# add this assay to the previously created Seurat object
Fixed_cells.i[["LOG2"]] <- log2_assay

# Validate that the object now contains multiple assays
Assays(Fixed_cells.i)
#Check the active assay 
DefaultAssay(Fixed_cells.i)
# Switch the default to LOG2
DefaultAssay(Fixed_cells.i) <- "LOG2"
DefaultAssay(Fixed_cells.i)


#Granule mature 
FeaturePlot(Fixed_cells.i, features= c("Prox1","Plk5", "Ntng1", "Rbfox3", "Calb1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Prox1","Plk5", "Ntng1", "Rbfox3", "Calb1"), slot ="counts")
#Granule immature
FeaturePlot(Fixed_cells.i, features= c("Ncam1","Prox1","Dcx", "Fxyd7","Rbfox3", "Calb2" ), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Ncam1","Prox1","Dcx", "Fxyd7","Rbfox3", "Calb2" ), slot ="counts")
#Neuroblast 2
FeaturePlot(Fixed_cells.i, features= c("Prox1","Dcx","Gal", "Sox11","Neurod1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Prox1","Dcx","Gal", "Sox11","Neurod1"), slot ="counts")
#Neuroblast1
FeaturePlot(Fixed_cells.i, features= c("Prox1","Dcx","Eomes", "Igfbpl1", "Calb2", "Neurog2", "Sox11", "Neurod4", "Tac2", "Neurod1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Prox1","Dcx","Eomes", "Igfbpl1", "Calb2", "Neurog2", "Sox11", "Neurod4", "Tac2", "Neurod1"), slot ="counts")
#nIPC 
FeaturePlot(Fixed_cells.i, features= c("Prox1","Dcx","Ascl", "Cdk1", "Tfap2c", "Eomes", "Igfbpl1","Neurog2", "Top2a", "Neurod4", "Aurkb", "Neurod1", "Ccnd2"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Prox1","Dcx","Ascl", "Cdk1", "Tfap2c", "Eomes", "Igfbpl1","Neurog2", "Top2a", "Neurod4", "Aurkb", "Neurod1", "Ccnd2"), slot ="counts")
#RGL
FeaturePlot(Fixed_cells.i, features= c("Tfap2c", "Sox9", "Lpar1", "Ascl1","Sox2", "Nes", "Vnn1","Hes5", "Rhcg", "Scl1a3", "Wnt8b", "Aldoc", "Apoe", "Id4", "Hopx", "GFAP"), slot ="counts")
VlnPlot(Fixed_cells.i, features= c("Tfap2c", "Sox9", "Lpar1", "Ascl1","Sox2", "Nes", "Vnn1","Hes5", "Rhcg", "Scl1a3", "Wnt8b", "Aldoc", "Apoe", "Id4", "Hopx", "GFAP"), slot ="counts")
#Astrocytes
FeaturePlot(Fixed_cells.i, features= c("S100b", "Fzd2", "Sox9","Sox2","Hes5","Scl1a3"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("S100b", "Fzd2", "Sox9","Sox2","Hes5","Scl1a3"), slot ="counts")
#OPC 
FeaturePlot(Fixed_cells.i, features= c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"), slot ="counts")
#Oligo 
FeaturePlot(Fixed_cells.i, features= c("Lpar1","Mbp", "Plp1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Lpar1","Mbp", "Plp1"), slot ="counts")
#Endothelial cells
FeaturePlot(Fixed_cells.i, features= c("Nes", "Prom1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Nes", "Prom1"), slot ="counts")
#Microglia
FeaturePlot(Fixed_cells.i, features= c("Csf1r" , "Cx3cr1"), slot ="counts",pt.size=1.5)
VlnPlot(Fixed_cells.i, features= c("Csf1r" , "Cx3cr1"), slot ="counts")

#------------------------------------------------Cluster investigation-----------------------------------------------------

#Playing with resolution
pilot.integrated <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=10)
pilot.integrated<- FindClusters(pilot.integrated, resolution = seq(from = 0.1, to = 1, by = 0.1), algorithm=1)
clustree(pilot.integrated, prefix = "integrated_snn_res.")
#Playing with K
pilot.integrated_K5 <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=5)
pilot.integrated_K7 <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=7)
pilot.integrated_K10 <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=10)
pilot.integrated_K15 <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=15)
pilot.integrated_K20 <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=20)

pilot.integrated_K5<- FindClusters(pilot.integrated_K5, resolution = 0.5, algorithm=1)
pilot.integrated_K7<- FindClusters(pilot.integrated_K7, resolution = 0.5, algorithm=1)
pilot.integrated_K10<- FindClusters(pilot.integrated_K10, resolution = 0.5, algorithm=1)
pilot.integrated_K15<- FindClusters(pilot.integrated_K15, resolution = 0.5, algorithm=1)
pilot.integrated_K20<- FindClusters(pilot.integrated_K20, resolution = 0.5, algorithm=1)



pilot.integrated[["K5"]] <- pilot.integrated_K5@meta.data[["seurat_clusters"]]
pilot.integrated[["K7"]] <- pilot.integrated_K7@meta.data[["seurat_clusters"]]
pilot.integrated[["K10"]] <- pilot.integrated_K10@meta.data[["seurat_clusters"]]
pilot.integrated[["K15"]] <- pilot.integrated_K15@meta.data[["seurat_clusters"]]
pilot.integrated[["K20"]] <- pilot.integrated_K20@meta.data[["seurat_clusters"]]

clustree(pilot.integrated, prefix = "K")

