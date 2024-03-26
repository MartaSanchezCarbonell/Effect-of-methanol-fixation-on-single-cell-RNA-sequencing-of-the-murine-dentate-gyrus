#Setup the Seurat Object
#Install pakages (only once)
install.packages('Seurat')
install.packages("devtools")
install.packages ("tidyverse")
install.packages("umap-learn")
install.packages("xlsx")
install.packages("rlang")
install.packages("ggpubr")
install.packages("car")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dittoSeq")
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
install.packages("Hmisc")
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
library(rlang)
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
# Load the pilot dataset
pilot.data<-readRDS("E:/Getting started with scRNAseq analysis/Pilot results/scd-22aeac/diagnostics/2021-10-11-JEN-AU.rds",refhook = NULL)

# In seurat we will not make use of the spike-ins, so remove them from pilot.data.rna before creating the Seurat object. 
## Calculate ERCC total reads and the abundance (%) on the raw counts before creating the Seurat object
ercc <- grep("ERCC-",rownames(pilot.data))
total_ercc_reads<-Matrix::colSums(pilot.data[ercc,])
percent.ERCC <- (Matrix::colSums(pilot.data[ercc,])/Matrix::colSums(pilot.data))*100

## Calculate mt total reads on the raw counts before creating a Seurat object
mito <- grep("^mt",rownames(pilot.data))
total_mito_reads<-Matrix::colSums(pilot.data[mito,])

## Convert the dataframe in a matrix and add Var_expression_per_cell, and Mean_expression_per_cell to the dataframe.
pilot.data.matrix<-as.matrix(pilot.data)
Var_expression_per_cell<-colVars(pilot.data.matrix)

Mean_expression_per_cell<-colMeans(pilot.data)

# When we create the seurat object we need to ask to remove the transcripts belonging to the ercc matrix and 
# add total_ercc_reads and percent.ERCC to object@meta.data in the total_ercc_reads and percent.ERCC column respectively
pilot <- CreateSeuratObject(counts = pilot.data[-ercc,], project = "Pilot analysis", min.cells = 0, min.features = 0)
pilot

# Remove chromosome lables
row.names(pilot.data)<-gsub("__.....","",row.names(pilot.data))
row.names(pilot.data)<-gsub("__....","",row.names(pilot.data))

pilot$orig.ident <- factor(pilot$orig.ident, levels=c("JEN-AU-s004","JEN-AU-s002"))

current.ids <- c("JEN-AU-s004","JEN-AU-s002")
new.ids <- c("Fresh_cells","Fixed_cells")
pilot@active.ident <- plyr::mapvalues(x = pilot@active.ident, from = current.ids, to = new.ids)
pilot$Condition <- Idents(pilot)

pilot$Condition <- factor(pilot$Condition, levels=c("Fresh_cells","Fixed_cells"))

# PRE-PROCESSING WORKFLOW

# Preparing all the extra columns in our surat Object
# QC and selecting cells for further analysis
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

# Histograms + cutoff lines
# Without y-axis fixed
hist(pilot$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fresh_cells$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fixed_cells$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5))
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

# WITH Y-AXIS FIXED (DONE BY PATRICIA) --> all graphs fixed according to pilot graphs.
hist(pilot$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5), ylim = c(0,10))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fresh_cells$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5), ylim = c(0,10))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")
hist(Fixed_cells$log10_nCount_RNA, breaks=500, xaxp=c(0,5,5), xlim = c(0,5), ylim = c(0,10))
abline(v=2.92942, col="blue")
abline(v=4.54407, col="blue")

hist(pilot$percent_mt, breaks=500, xlim = c(0,100), ylim = c(0,8))
hist(Fresh_cells$percent_mt,breaks=500, xlim = c(0,100), ylim = c(0,8))
hist(Fixed_cells$percent_mt,breaks=500, xlim = c(0,100), ylim = c(0,8))

hist(pilot$log10_nFeature_RNA, breaks=500, xlim = c(0,4), ylim = c(0,15))
abline(v=2.69897, col="blue")
hist(Fresh_cells$log10_nFeature_RNA,breaks=500, xlim = c(0,4), ylim = c(0,15))
abline(v=2.69897, col="blue")
hist(Fixed_cells$log10_nFeature_RNA,breaks=500, xlim = c(0,4), ylim = c(0,15))
abline(v=2.69897, col="blue")

hist(pilot$total_ercc_reads, breaks=500, xlim = c(0,2500), ylim = c(0,8))
abline(v=500, col="blue")
hist(Fresh_cells$total_ercc_reads,breaks=500, xlim = c(0,2500), ylim = c(0,8))
abline(v=500, col="blue")
hist(Fixed_cells$total_ercc_reads,breaks=500, xlim = c(0,2500), ylim = c(0,8))
abline(v=500, col="blue")

hist(pilot$ratio_mol_gen, breaks=500, xlim = c(0,55), ylim = c(0,50))
abline(v=1.2, col="blue")
hist(Fresh_cells$ratio_mol_gen,breaks=500, xlim = c(0,55), ylim = c(0,50))
abline(v=1.2, col="blue")
hist(Fixed_cells$ratio_mol_gen,breaks=500, xlim = c(0,55), ylim = c(0,50))
abline(v=1.2, col="blue")

# Filtering
Fresh_cells<-subset(Fresh_cells, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
Fixed_cells<-subset(Fixed_cells, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
pilot<-subset(pilot, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
Fresh_cells
Fixed_cells
pilot

# NORMALIZATION
Fresh_cells<-SCTransform(Fresh_cells,new.assay.name="SCT", variable.features.n=2000)
Fresh_cells

Fixed_cells<-SCTransform(Fixed_cells,new.assay.name="SCT", variable.features.n=2000)
Fixed_cells

#Integration --> DO NOT DO IT FOR THE NEXT PART (GENE NUMBER RETAINED AFTER EXPRESSION FILTERING), BUT DO IT FOR THE SCATTER PLOTS GRAPHS.
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

# Evaluate dropout frequency over the entire transcriptome
# Set a series (from less restrictive to high restrictive) of increasing gene expression level thresholds (nCount_RNA) 
# for defining detected genes -- BULK ANALYSIS

# RNA@counts

# Fresh_cells
DefaultAssay(Fresh_cells) <- "RNA"

pilot.bulk_RNA_counts_Fresh_cells<-AverageExpression(Fresh_cells, slot= "counts")$RNA

# I cannot eliminate the --chr fraction because there is no pilot.bulk_RNA_counts_Fresh_cells$data.
# I cannot neither change the identities because pilot.bulk_RNA_counts_Fresh_cells is not a Seurat object.

# I created a dataframe of the pilot because in the boxplot formula  data = dataframe --> DO IT JUST ONCE!
write.csv2(pilot.bulk_RNA_counts_Fresh_cells, (paste("E:/Getting started with scRNAseq analysis/Pilot results/7 Evaluate dropouts/pilot.bulk_RNA_counts_Fresh_cells.csv")), 
                                                         quote = FALSE, col.names = TRUE, row.names = TRUE)

# Create a dataframe with the RNA data (https://rdrr.io/r/utils/read.table.html)
pilot.bulk_RNA_counts_Fresh_cells <- as.data.frame(read.csv2(file = "E:/Getting started with scRNAseq analysis/Pilot results/7 Evaluate dropouts/pilot.bulk_RNA_counts_Fresh_cells.csv",
                                                        header = TRUE, dec = ",", row.names = 1))

boxplot(pilot.bulk_RNA_counts_Fresh_cells, data = pilot.bulk_RNA_counts_Fresh_cells, 
        xlab = "Fresh_cells", ylab = "Average expression", col = c("#F8766D"))

median(pilot.bulk_RNA_counts_Fresh_cells$all)
mean(pilot.bulk_RNA_counts_Fresh_cells$all)

#--------------------------------------------------------------------------------------

## SCT@data

# Fresh_cells
DefaultAssay(Fresh_cells) <- "SCT"

pilot.bulk_SCT_data_Fresh_cells<-AverageExpression(Fresh_cells, slot= "data")$SCT

pilot.bulk_SCT_data_Fresh_cells_df<-as.data.frame(pilot.bulk_SCT_data_Fresh_cells)

boxplot(pilot.bulk_SCT_data_Fresh_cells_df, data = pilot.bulk_SCT_data_Fresh_cells_df, 
        xlab = "Fresh_cells", ylab = "Average expression", col = c("#F8766D"))

max(pilot.bulk_SCT_data_Fresh_cells$all)
min(pilot.bulk_SCT_data_Fresh_cells$all)


#--------------------------------------------------------------------------------------

pilot.bulk_RNA_counts_Fresh_cells_dropouts_0.1<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 0.1)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_0.1, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_0.1, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_0.5<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 0.5)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_0.5, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_0.5, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_1<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 1)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_1, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_1, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_1.5<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 1.5)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_1.5, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_1.5, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_2<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 2)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_2, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_2, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_3<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 3)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_3, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_3, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_6<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 6)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_6, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_6, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_8<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 8)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_8, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_8, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_10<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 10)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_10, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_10, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_12.5<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 12.5)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_12.5, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_12.5, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_15<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 15)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_15, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_15, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_27<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 27)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_27, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_27, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_35<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 35)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_35, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_35, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

pilot.bulk_RNA_counts_Fresh_cells_dropouts_50<-subset.data.frame(pilot.bulk_RNA_counts_Fresh_cells, subset = all > 50)

boxplot(pilot.bulk_RNA_counts_Fresh_cells_dropouts_50, data = pilot.bulk_RNA_counts_Fresh_cells_dropouts_50, 
        xlab = "Condition", ylab = "Average expression", col = c("#F8766D"))

# Fixed_cells
DefaultAssay(Fixed_cells) <- "RNA"

pilot.bulk_RNA_counts_Fixed_cells<-AverageExpression(Fixed_cells, slot= "counts")$RNA

# I cannot eliminate the --chr fraction because there is no pilot.bulk_RNA_counts_Fixed_cells$data.
# I cannot neither change the identities because pilot.bulk_RNA_counts_Fixed_cells is not a Seurat object.

# I created a dataframe of the pilot because in the boxplot formula  data = dataframe --> DO IT JUST ONCE!
write.csv2(pilot.bulk_RNA_counts_Fixed_cells, (paste("E:/Getting started with scRNAseq analysis/Pilot results/7 Evaluate dropouts/pilot.bulk_RNA_counts_Fixed_cells.csv")), 
           quote = FALSE, col.names = TRUE, row.names = TRUE)

# Create a dataframe with the RNA data (https://rdrr.io/r/utils/read.table.html)
pilot.bulk_RNA_counts_Fixed_cells <- as.data.frame(read.csv2(file = "E:/Getting started with scRNAseq analysis/Pilot results/7 Evaluate dropouts/pilot.bulk_RNA_counts_Fixed_cells.csv",
                                                             header = TRUE, dec = ",", row.names = 1))

boxplot(pilot.bulk_RNA_counts_Fixed_cells, data = pilot.bulk_RNA_counts_Fixed_cells, xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

median(pilot.bulk_RNA_counts_Fixed_cells$all)
mean(pilot.bulk_RNA_counts_Fixed_cells$all)

#--------------------------------------------------------------------------------------

## SCT@data

# Fixed_cells
DefaultAssay(Fixed_cells) <- "SCT"

pilot.bulk_SCT_data_Fixed_cells<-AverageExpression(Fixed_cells, slot= "data")$SCT

pilot.bulk_SCT_data_Fixed_cells_df<-as.data.frame(pilot.bulk_SCT_data_Fixed_cells)

boxplot(pilot.bulk_SCT_data_Fixed_cells_df, data = pilot.bulk_SCT_data_Fixed_cells_df, 
        xlab = "Fixed_cells", ylab = "Average expression", col = c("#00B8E7"))

max(pilot.bulk_SCT_data_Fresh_cells$all)
min(pilot.bulk_SCT_data_Fresh_cells$all)


#--------------------------------------------------------------------------------------

pilot.bulk_RNA_counts_Fixed_cells_dropouts_0.1<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 0.1)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_0.1, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_0.1, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_0.5<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 0.5)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_0.5, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_0.5, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_1<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 1)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_1, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_1, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_1.5<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 1.5)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_1.5, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_1.5, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))


pilot.bulk_RNA_counts_Fixed_cells_dropouts_2<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 2)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_2, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_2, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_3<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 3)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_3, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_3, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_6<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 6)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_6, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_6, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_8<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 8)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_8, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_8, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_10<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 10)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_10, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_10, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_12.5<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 12.5)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_12.5, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_12.5, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_15<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 15)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_15, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_15, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_27<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 27)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_27, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_27, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_35<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 35)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_35, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_35, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))

pilot.bulk_RNA_counts_Fixed_cells_dropouts_50<-subset.data.frame(pilot.bulk_RNA_counts_Fixed_cells, subset = all > 50)

boxplot(pilot.bulk_RNA_counts_Fixed_cells_dropouts_50, data = pilot.bulk_RNA_counts_Fixed_cells_dropouts_50, 
        xlab = "Condition", ylab = "Average expression", col = c("#00B8E7"))



DefaultAssay(pilot.integrated) <- "RNA"

# log2(+1)=log2_1

pilot.bulk_RNA_counts_log2_1<-log2(1+(AverageExpression(pilot.integrated, group.by="Condition", slot= "counts"))$RNA)

# I created a dataframe of the pilot because in the boxplot formula  data = dataframe --> DO IT JUST ONCE!
write.csv2(pilot.bulk_RNA_counts_log2_1, (paste("E:/Getting started with scRNAseq analysis/Pilot results/7 Evaluate dropouts/pilot.bulk_RNA_counts_log2_1.csv")), 
           quote = FALSE, col.names = TRUE, row.names = TRUE)

# Create a dataframe with the RNA data (https://rdrr.io/r/utils/read.table.html)
pilot.bulk_RNA_counts_log2_1 <- as.data.frame(read.csv2(file = "E:/Getting started with scRNAseq analysis/Pilot results/7 Evaluate dropouts/pilot.bulk_RNA_counts_log2_1.csv",
                                                             header = TRUE, dec = ",", row.names = 1))

# Low gene expression
# pilot.bulk_RNA_counts_log2_1_dropouts_0.1<-subset.data.frame(pilot.bulk_RNA_counts_log2_1, subset = Fresh_cells < 0.13750 & Fixed_cells < 0.13750)
pilot.bulk_RNA_counts_log2_1_dropouts_1<-subset.data.frame(pilot.bulk_RNA_counts_log2_1, subset = Fresh_cells < 1 & Fixed_cells < 1)

# Scatter plot

pilot.bulk_RNA_counts_log2_1_dropouts_1<-data.frame(pilot.bulk_RNA_counts_log2_1_dropouts_1)

p1 <- ggplot(pilot.bulk_RNA_counts_log2_1_dropouts_1, aes(Fixed_cells, Fresh_cells)) + geom_point(size=1) +
  theme_light() +geom_abline(intercept = 0, slope = 1, color = "blue")+ 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("log2(average expression+1) in Fixed_cells")+ xlim(0,3)+
  theme(axis.title.x = element_text(family="Calibri Light" , size=10))+
  ylab("log2(average expression+1) in Fresh_cells")+ ylim(0,3)+
  theme(axis.title.y = element_text(family="Calibri Light" , size=10))

p1

mean(pilot.bulk_RNA_counts_log2_1_dropouts_1$Fresh_cells)
mean(pilot.bulk_RNA_counts_log2_1_dropouts_1$Fixed_cells)

# High gene expression
pilot.bulk_RNA_counts_log2_1_dropouts_8<-subset.data.frame(pilot.bulk_RNA_counts_log2_1, subset = Fresh_cells > 3.16993 & Fixed_cells > 3.16993)

# Scatter plot

pilot.bulk_RNA_counts_log2_1_dropouts_8<-data.frame(pilot.bulk_RNA_counts_log2_1_dropouts_8)

p2 <- ggplot(pilot.bulk_RNA_counts_log2_1_dropouts_8, aes(Fixed_cells, Fresh_cells)) + geom_point(size=1) +
  theme_light() +geom_abline(intercept = 0, slope = 1, color = "blue")+ 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("log2(average expression+1) in Fixed_cells")+xlim(2,10)+
  theme(axis.title.x = element_text(family="Calibri Light" , size=10))+
  ylab("log2(average expression+1) in Fresh_cells")+ylim(2,10)+
  theme(axis.title.y = element_text(family="Calibri Light" , size=10))

p2

mean(pilot.bulk_RNA_counts_log2_1_dropouts_8$Fresh_cells)
mean(pilot.bulk_RNA_counts_log2_1_dropouts_8$Fixed_cells)

# No gene expression filtering
# Scatter plot

pilot.bulk_RNA_counts_log2_1<-data.frame(pilot.bulk_RNA_counts_log2_1)

p3 <- ggplot(pilot.bulk_RNA_counts_log2_1, aes(Fixed_cells, Fresh_cells)) + geom_point(size=1) +
  theme_light() +geom_abline(intercept = 0, slope = 1, color = "blue")+ 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("log2(average expression+1) in Fixed_cells")+ xlim(0,8)+
  theme(axis.title.x = element_text(family="Calibri Light" , size=10))+
  ylab("log2(average expression+1) in Fresh_cells")+ ylim(0,8)+
  theme(axis.title.y = element_text(family="Calibri Light" , size=10))

p3

#SCT@data

DefaultAssay(pilot.integrated) <- "SCT"

#log2(+1)=log2_1

pilot.bulk_SCT_data_log2_1<-log2(1+(AverageExpression(pilot.integrated, group.by="Condition", slot= "data"))$SCT)

# I created a dataframe of the pilot because in the boxplot formula  data = datafram --> DO IT JUST ONCE!
write.csv2(pilot.bulk_SCT_data_log2_1, (paste("E:/Getting started with scRNAseq analysis/Pilot results/7 Evaluate dropouts/pilot.bulk_SCT_data_log2_1.csv")), 
           quote = FALSE, col.names = TRUE, row.names = TRUE)

# Create a dataframe with the RNA data (https://rdrr.io/r/utils/read.table.html)
pilot.bulk_SCT_data_log2_1 <- as.data.frame(read.csv2(file = "E:/Getting started with scRNAseq analysis/Pilot results/7 Evaluate dropouts/pilot.bulk_SCT_data_log2_1.csv",
                                                        header = TRUE, dec = ",", row.names = 1))

# Low gene expression
pilot.bulk_SCT_data_log2_1_dropouts_1<-subset.data.frame(pilot.bulk_SCT_data_log2_1, subset = Fresh_cells < 1 & Fixed_cells < 1)

# Scatter plot

pilot.bulk_SCT_data_log2_1_dropouts_1<-data.frame(pilot.bulk_SCT_data_log2_1_dropouts_1)

p3 <- ggplot(pilot.bulk_SCT_data_log2_1_dropouts_1, aes(Fixed_cells, Fresh_cells)) + geom_point(size=1) + 
  theme_light() +geom_abline(intercept = 0, slope = 1, color = "blue")+ 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("log2(average expression+1) in Fixed_cells")+xlim(0,3)+
  theme(axis.title.x = element_text(family="Calibri Light" , size=10))+
  ylab("log2(average expression+1) in Fresh_cells")+ylim(0,3)+
  theme(axis.title.y = element_text(family="Calibri Light" , size=10))

p3

mean(pilot.bulk_SCT_data_log2_1_dropouts_1$Fresh_cells)
mean(pilot.bulk_SCT_data_log2_1_dropouts_1$Fixed_cells)

# High gene expression
# Low gene expression
pilot.bulk_SCT_data_log2_1_dropouts_8<-subset.data.frame(pilot.bulk_SCT_data_log2_1, subset = Fresh_cells > 3.16993 & Fixed_cells > 3.16993)

# Scatter plot

pilot.bulk_SCT_data_log2_1_dropouts_8<-data.frame(pilot.bulk_SCT_data_log2_1_dropouts_8)

p4 <- ggplot(pilot.bulk_SCT_data_log2_1_dropouts_8, aes(Fixed_cells, Fresh_cells)) + geom_point(size=1) + 
  theme_light() +geom_abline(intercept = 0, slope = 1, color = "blue")+
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("log2(average expression+1) in Fixed_cells")+xlim(2,10)+
  theme(axis.title.x = element_text(family="Calibri Light" , size=10))+
  ylab("log2(average expression+1) in Fresh_cells")+ylim(2,10)+
  theme(axis.title.y = element_text(family="Calibri Light" , size=10))

p4

mean(pilot.bulk_SCT_data_log2_1_dropouts_8$Fresh_cells)
mean(pilot.bulk_SCT_data_log2_1_dropouts_8$Fixed_cells)

# No gene expression filtering
# Scatter plot

pilot.bulk_SCT_data_log2_1<-data.frame(pilot.bulk_SCT_data_log2_1)

p5 <- ggplot(pilot.bulk_SCT_data_log2_1, aes(Fixed_cells, Fresh_cells)) + geom_point(size=1) +
  theme_light() +geom_abline(intercept = 0, slope = 1, color = "blue")+ 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("log2(average expression+1) in Fixed_cells")+ xlim(0,8)+
  theme(axis.title.x = element_text(family="Calibri Light" , size=10))+
  ylab("log2(average expression+1) in Fresh_cells")+ ylim(0,8)+
  theme(axis.title.y = element_text(family="Calibri Light" , size=10))

p5



