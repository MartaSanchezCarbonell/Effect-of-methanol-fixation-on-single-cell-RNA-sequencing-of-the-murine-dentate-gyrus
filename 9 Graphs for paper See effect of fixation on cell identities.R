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


#After cluster annotation with Linnarsson's database
#Assigning cell type identity to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("Astrocytes +  RGL", "nIPCs + Neuroblasts", "Immature neurons", 
                     "Mature neurons","Microglia","Oligodendrocytes", "OPCs", "Endothelial")
pilot.integrated@active.ident <- plyr::mapvalues(x = pilot.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
pilot.integrated$CellType <- Idents(pilot.integrated)

#Save
saveRDS(pilot.integrated, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/9 See effect of fixation on cell identities/R objects/Graphs_for_paper_pilot.integrated.rds")
#Load
pilot.integrated<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/9 See effect of fixation on cell identities/R objects/Graphs_for_paper_pilot.integrated.rds")


table(pilot.integrated@meta.data[["CellType"]],pilot.integrated@meta.data[["Condition"]])


#UMAP
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15, n.neighbors = 100)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pilot.integrated, reduction = "umap", label=FALSE, pt.size=2.3)



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
mdata <- as.data.frame(pilot.integrated@meta.data)
Cell_type<-mdata$CellType
mdata2<-as.data.frame(Cell_type)
mdata2$Condition<-mdata$Condition

Fresh_df<-subset(mdata2, Condition=="Fresh_cells")
Fixed_df<-subset(mdata2, Condition=="Fixed_cells")
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

ggarrange(Fr,Fi)



time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fresh_per_cent, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/9 See effect of fixation on cell identities/", paste(time.name), "Percetages_CellType_Condition.xlsx"))
           , sheetName = "Fresh", col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(Fixed_per_cent, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/9 See effect of fixation on cell identities/", paste(time.name), "Percetages_CellType_Condition.xlsx"))
           , sheetName = "Fixed", col.names = TRUE, row.names = TRUE, append = T)


#---------------------------top 10 differentially expressed genes heat map with cells labelled with cluster and condition--------------------------------------------------------
#Load
pilot.integrated<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/9 See effect of fixation on cell identities/R objects/Graphs_for_paper_pilot.integrated.rds")

#Lists of differential expressed genes (clusters)
pilot.markers <- FindAllMarkers(object=pilot.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.001, test.use= "wilcox",assay="RNA")


#Differential expressed genes heatmap (clusters)
pilot.markers %>%
  group_by(cluster) %>%
  slice_max(n =10, order_by = avg_log2FC) -> top10_log2FC
top10_log2FC_HM<-DoHeatmap(pilot.integrated, assay="RNA", features = top10_log2FC$gene, label=FALSE)
top10_log2FC_HM + labs(title="top10_log2FC
                       ")

#--------Heatmap-------------
#https://github.com/elliefewings/DoMultiBarHeatmap
#devtools::install_github("elliefewings/DoMultiBarHeatmap")

#Load library
library(DoMultiBarHeatmap)
#Cell order
Clusters<-(as.data.frame(pilot.integrated$CellType))
colnames(Clusters) <- c('CellType')
Clusters$Cells<-rownames(Clusters)
Clusters$Num<-pilot.integrated$seurat_clusters
Clusters<-Clusters[order(Clusters$Num,decreasing=F),]
Cell_order<-Clusters$Cells

top10_log2FC_HM<-DoMultiBarHeatmap(pilot.integrated, assay = 'SCT', cells=factor(Cell_order), features = top10_log2FC$gene, group.by='Condition', additional.group.by = 'CellType')

top10_log2FC_HM + labs(title="top10_log2FC")

top10_log2FC_HM_2<-DoMultiBarHeatmap(pilot.integrated, assay = 'integrated', features = top10_log2FC$gene, group.by='seurat_clusters', additional.group.by = 'Condition')

top10_log2FC_HM_2 + labs(title="top10_log2FC_centered_and_corrected_Pearson_residuals")
