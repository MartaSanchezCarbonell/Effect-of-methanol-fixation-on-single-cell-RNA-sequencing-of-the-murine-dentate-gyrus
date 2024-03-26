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








#RNA@counts




DefaultAssay(pilot.integrated) <- "RNA"

      #log1p

pilot.bulk_RNA_counts_log1p<-log1p(AverageExpression(pilot.integrated, group.by="Condition", slot= "counts")$RNA)
pilot.bulk_RNA_counts_log1p<-data.frame(pilot.bulk_RNA_counts_log1p)

            #Plot density distribution

ggdensity(pilot.bulk_RNA_counts_log1p$Fresh_cells, 
          main = "Pilot.bulk_RNA_counts$Fresh_cells",
          xlab = "Counts")
ggdensity(pilot.bulk_RNA_counts_log1p$Fixed_cells, 
          main = "Pilot.bulk_RNA_counts$Fixed_cells",
          xlab = "Counts")

            #Normal quantile plot

qqnorm(pilot.bulk_RNA_counts_log1p$Fresh_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_RNA_counts_log1p$Fresh_cells quantiles")
qqline(pilot.bulk_RNA_counts_log1p$Fresh_cells, col = "steelblue", lwd = 2)
qqnorm(pilot.bulk_RNA_counts_log1p$Fixed_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_RNA_counts_log1p$Fixed_cells quantiles")
qqline(pilot.bulk_RNA_counts_log1p$Fixed_cells, col = "steelblue", lwd = 2)

            #Plot scatter plot

p1 <- ggplot(pilot.bulk_RNA_counts_log1p, aes(Fixed_cells,Fresh_cells)) + geom_point(size=2) + ggtitle("Comparation")+ 
              theme_light() +geom_smooth(method=lm, se =TRUE) + 
              theme(plot.title = element_text(family="Calibri Light" , size=20))+
              xlab("Fixed cells (ln(1+counts)")+
              theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
              ylab("Fresh cells (ln(1+counts)")+
              theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p1
            #Run Spearman’s correlation test


cor.test(pilot.bulk_RNA_counts_log1p$Fresh_cells, pilot.bulk_RNA_counts_log1p$Fixed_cells, method=c("spearman"))
cor.test(pilot.bulk_RNA_counts_log1p$Fresh_cells, pilot.bulk_RNA_counts_log1p$Fixed_cells, method=c("pearson"))

      #log2

pilot.bulk_RNA_counts_log2<-log2(AverageExpression(pilot.integrated, group.by="Condition", slot= "counts")$RNA)
pilot.bulk_RNA_counts_log2<-data.frame(pilot.bulk_RNA_counts_log2)

            #Plot density distribution

ggdensity(pilot.bulk_RNA_counts_log2$Fresh_cells, 
          main = "Pilot.bulk_RNA_counts$Fresh_cells",
          xlab = "Counts")
ggdensity(pilot.bulk_RNA_counts_log2$Fixed_cells, 
          main = "Pilot.bulk_RNA_counts$Fixed_cells",
          xlab = "Counts")

            #Normal quantile plot

qqnorm(pilot.bulk_RNA_counts_log2$Fresh_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_RNA_counts_log2$Fresh_cells quantiles",na.rm=TRUE)
qqline(pilot.bulk_RNA_counts_log2$Fresh_cells, col = "steelblue", lwd = 2)
qqnorm(pilot.bulk_RNA_counts_log2$Fixed_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_RNA_counts_log2$Fixed_cells quantiles")
qqline(pilot.bulk_RNA_counts_log2$Fixed_cells, col = "steelblue", lwd = 2)

            #Plot scatter plot

p1 <- ggplot(pilot.bulk_RNA_counts_log2, aes(Fixed_cells,Fresh_cells)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() +geom_smooth(method=lm, se =TRUE) + 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (log2(counts)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (log2(counts)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p1
            #Run Spearman’s correlation test


cor.test(pilot.bulk_RNA_counts_log2$Fresh_cells, pilot.bulk_RNA_counts_log2$Fixed_cells, method=c("spearman"))
cor.test(pilot.bulk_RNA_counts_log2$Fresh_cells, pilot.bulk_RNA_counts_log2$Fixed_cells, method=c("pearson"))


      #log2(+1)=log2_1

pilot.bulk_RNA_counts_log2_1<-log2(1+(AverageExpression(pilot.integrated, group.by="Condition", slot= "counts"))$RNA)
pilot.bulk_RNA_counts_log2_1<-data.frame(pilot.bulk_RNA_counts_log2_1)

            #Plot density distribution

ggdensity(pilot.bulk_RNA_counts_log2_1$Fresh_cells, 
          main = "Pilot.bulk_RNA_counts$Fresh_cells",
          xlab = "Counts")
ggdensity(pilot.bulk_RNA_counts_log2_1$Fixed_cells, 
          main = "Pilot.bulk_RNA_counts$Fixed_cells",
          xlab = "Counts")

            #Normal quantile plot

qqnorm(pilot.bulk_RNA_counts_log2_1$Fresh_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_RNA_counts_log2_1$Fresh_cells quantiles",na.rm=TRUE)
qqline(pilot.bulk_RNA_counts_log2_1$Fresh_cells, col = "steelblue", lwd = 2)
qqnorm(pilot.bulk_RNA_counts_log2_1$Fixed_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_RNA_counts_log2_1$Fixed_cells quantiles")
qqline(pilot.bulk_RNA_counts_log2_1$Fixed_cells, col = "steelblue", lwd = 2)

            #Plot scatter plot

p1 <- ggplot(pilot.bulk_RNA_counts_log2_1, aes(Fixed_cells,Fresh_cells)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() +geom_smooth(method=lm, se =TRUE) + 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (log2(counts+1)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (log2(counts+1)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p1

p2 <- ggplot(pilot.bulk_RNA_counts_log2_1, aes(Fixed_cells,Fresh_cells)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() + geom_abline(intercept = 0, slope = 1, color = "blue")+ xlim(0, 9)+ylim(0, 9)+
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (log2(counts+1)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (log2(counts+1)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p2



            #Run Spearman’s correlation test


cor.test(pilot.bulk_RNA_counts_log2_1$Fresh_cells, pilot.bulk_RNA_counts_log2_1$Fixed_cells, method=c("spearman"))
cor.test(pilot.bulk_RNA_counts_log2_1$Fresh_cells, pilot.bulk_RNA_counts_log2_1$Fixed_cells, method=c("pearson"))


      #log10(+1)=log10_1

pilot.bulk_RNA_counts_log10_1<-log10(1+(AverageExpression(pilot.integrated, group.by="Condition", slot= "counts"))$RNA)
pilot.bulk_RNA_counts_log10_1<-data.frame(pilot.bulk_RNA_counts_log10_1)

            #Plot density distribution

ggdensity(pilot.bulk_RNA_counts_log10_1$Fresh_cells, 
          main = "Pilot.bulk_RNA_counts$Fresh_cells",
          xlab = "Counts")
ggdensity(pilot.bulk_RNA_counts_log10_1$Fixed_cells, 
          main = "Pilot.bulk_RNA_counts$Fixed_cells",
          xlab = "Counts")

            #Normal quantile plot

qqnorm(pilot.bulk_RNA_counts_log10_1$Fresh_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_RNA_counts_log10_1$Fresh_cells quantiles",na.rm=TRUE)
qqline(pilot.bulk_RNA_counts_log10_1$Fresh_cells, col = "steelblue", lwd = 2)
qqnorm(pilot.bulk_RNA_counts_log10_1$Fixed_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_RNA_counts_log10_1$Fixed_cells quantiles")
qqline(pilot.bulk_RNA_counts_log10_1$Fixed_cells, col = "steelblue", lwd = 2)

            #Plot scatter plot

p1 <- ggplot(pilot.bulk_RNA_counts_log10_1, aes(Fixed_cells,Fresh_cells)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() +geom_smooth(method=lm, se =TRUE) + 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (log10(counts+1)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (log10(counts+1)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p1
              #Run Spearman’s correlation test


cor.test(pilot.bulk_RNA_counts_log10_1$Fresh_cells, pilot.bulk_RNA_counts_log10_1$Fixed_cells, method=c("spearman"))
cor.test(pilot.bulk_RNA_counts_log10_1$Fresh_cells, pilot.bulk_RNA_counts_log10_1$Fixed_cells, method=c("pearson"))








#SCT@data







DefaultAssay(pilot.integrated) <- "SCT"

      #log1p

pilot.bulk_SCT_data_log1p<-log1p(AverageExpression(pilot.integrated, group.by="Condition", slot= "data")$SCT)
pilot.bulk_SCT_data_log1p<-data.frame(pilot.bulk_SCT_data_log1p)

            #Plot density distribution

ggdensity(pilot.bulk_SCT_data_log1p$Fresh_cells, 
          main = "Pilot.bulk_SCT_data$Fresh_cells",
          xlab = "Counts")
ggdensity(pilot.bulk_SCT_data_log1p$Fixed_cells, 
          main = "Pilot.bulk_SCT_data$Fixed_cells",
          xlab = "Counts")

            #Normal quantile plot

qqnorm(pilot.bulk_SCT_data_log1p$Fresh_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_SCT_data_log1p$Fresh_cells quantiles")
qqline(pilot.bulk_SCT_data_log1p$Fresh_cells, col = "steelblue", lwd = 2)
qqnorm(pilot.bulk_SCT_data_log1p$Fixed_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_SCT_data_log1p$Fixed_cells quantiles")
qqline(pilot.bulk_SCT_data_log1p$Fixed_cells, col = "steelblue", lwd = 2)

#Plot scatter plot

p1 <- ggplot(pilot.bulk_SCT_data_log1p, aes(Fixed_cells,Fresh_cells)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() +geom_smooth(method=lm, se =TRUE) + 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (ln(1+counts)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (ln(1+counts)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p1

              #Run Spearman’s correlation test


cor.test(pilot.bulk_SCT_data_log1p$Fresh_cells, pilot.bulk_SCT_data_log1p$Fixed_cells, method=c("spearman"))
cor.test(pilot.bulk_SCT_data_log1p$Fresh_cells, pilot.bulk_SCT_data_log1p$Fixed_cells, method=c("pearson"))

      #log2

pilot.bulk_SCT_data_log2<-log2(AverageExpression(pilot.integrated, group.by="Condition", slot= "data")$SCT)
pilot.bulk_SCT_data_log2<-data.frame(pilot.bulk_SCT_data_log2)

            #Plot density distribution

ggdensity(pilot.bulk_SCT_data_log2$Fresh_cells, 
          main = "Pilot.bulk_SCT_data$Fresh_cells",
          xlab = "Counts")
ggdensity(pilot.bulk_SCT_data_log2$Fixed_cells, 
          main = "Pilot.bulk_SCT_data$Fixed_cells",
          xlab = "Counts")

            #Normal quantile plot

qqnorm(pilot.bulk_SCT_data_log2$Fresh_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_SCT_data_log2$Fresh_cells quantiles")
qqline(pilot.bulk_SCT_data_log2$Fresh_cells, col = "steelblue", lwd = 2)
qqnorm(pilot.bulk_SCT_data_log2$Fixed_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_SCT_data_log2$Fixed_cells quantiles")
qqline(pilot.bulk_SCT_data_log2$Fixed_cells, col = "steelblue", lwd = 2)

            #Plot scatter plot

p1 <- ggplot(pilot.bulk_SCT_data_log2, aes(Fixed_cells,Fresh_cells)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() +geom_smooth(method=lm, se =TRUE) + 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (log2(counts)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (log2(counts)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p1

            #Run Spearman’s correlation test


cor.test(pilot.bulk_SCT_data_log2$Fresh_cells, pilot.bulk_SCT_data_log2$Fixed_cells, method=c("spearman"))
cor.test(pilot.bulk_SCT_data_log2$Fresh_cells, pilot.bulk_SCT_data_log2$Fixed_cells, method=c("pearson"))




      #log2(+1)=log2_1

pilot.bulk_SCT_data_log2_1<-log2(1+(AverageExpression(pilot.integrated, group.by="Condition", slot= "data"))$SCT)
pilot.bulk_SCT_data_log2_1<-data.frame(pilot.bulk_SCT_data_log2_1)

            #Plot density distribution

ggdensity(pilot.bulk_SCT_data_log2_1$Fresh_cells, 
          main = "Pilot.bulk_SCT_data$Fresh_cells",
          xlab = "Counts")
ggdensity(pilot.bulk_SCT_data_log2_1$Fixed_cells, 
          main = "Pilot.bulk_SCT_data$Fixed_cells",
          xlab = "Counts")

            #Normal quantile plot

qqnorm(pilot.bulk_SCT_data_log2_1$Fresh_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_SCT_data_log2_1$Fresh_cells quantiles",na.rm=TRUE)
qqline(pilot.bulk_SCT_data_log2_1$Fresh_cells, col = "steelblue", lwd = 2)
qqnorm(pilot.bulk_SCT_data_log2_1$Fixed_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_SCT_data_log2_1$Fixed_cells quantiles")
qqline(pilot.bulk_SCT_data_log2_1$Fixed_cells, col = "steelblue", lwd = 2)

            #Plot scatter plot

p1 <- ggplot(pilot.bulk_SCT_data_log2_1, aes(Fixed_cells,Fresh_cells)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() +geom_smooth(method=lm, se =TRUE) + 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (log2(counts+1)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (log2(counts+1)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p1

            #Run Spearman’s correlation test


cor.test(pilot.bulk_SCT_data_log2_1$Fresh_cells, pilot.bulk_SCT_data_log2_1$Fixed_cells, method=c("spearman"))
cor.test(pilot.bulk_SCT_data_log2_1$Fresh_cells, pilot.bulk_SCT_data_log2_1$Fixed_cells, method=c("pearson"))




      #log10(+1)=log10_1

pilot.bulk_SCT_data_log10_1<-log10(1+(AverageExpression(pilot.integrated, group.by="Condition", slot= "data"))$SCT)
pilot.bulk_SCT_data_log10_1<-data.frame(pilot.bulk_SCT_data_log10_1)

            #Plot density distribution

ggdensity(pilot.bulk_SCT_data_log10_1$Fresh_cells, 
          main = "Pilot.bulk_SCT_data$Fresh_cells",
          xlab = "Counts")
ggdensity(pilot.bulk_SCT_data_log10_1$Fixed_cells, 
          main = "Pilot.bulk_SCT_data$Fixed_cells",
          xlab = "Counts")

            #Normal quantile plot

qqnorm(pilot.bulk_SCT_data_log10_1$Fresh_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_SCT_data_log10_1$Fresh_cells quantiles",na.rm=TRUE)
qqline(pilot.bulk_SCT_data_log10_1$Fresh_cells, col = "steelblue", lwd = 2)
qqnorm(pilot.bulk_SCT_data_log10_1$Fixed_cells, pch = 1, frame = TRUE, xlab = "Normal distribution quantiles", ylab = "pilot.bulk_SCT_data_log10_1$Fixed_cells quantiles")
qqline(pilot.bulk_SCT_data_log10_1$Fixed_cells, col = "steelblue", lwd = 2)

            #Plot scatter plot

p1 <- ggplot(pilot.bulk_SCT_data_log10_1, aes(Fixed_cells,Fresh_cells)) + geom_point(size=2) + ggtitle("Comparation")+ 
  theme_light() +geom_smooth(method=lm, se =TRUE) + 
  theme(plot.title = element_text(family="Calibri Light" , size=20))+
  xlab("Fixed cells (log10(counts+1)")+
  theme(axis.title.x = element_text(family="Calibri Light" , size=20))+
  ylab("Fresh cells (log10(counts+1)")+
  theme(axis.title.y = element_text(family="Calibri Light" , size=20))

p1

            #Run Spearman’s correlation test


cor.test(pilot.bulk_SCT_data_log10_1$Fresh_cells, pilot.bulk_SCT_data_log10_1$Fixed_cells, method=c("spearman"))
cor.test(pilot.bulk_SCT_data_log10_1$Fresh_cells, pilot.bulk_SCT_data_log10_1$Fixed_cells, method=c("pearson"))













#To see the impact of preservation on cellular and pseudo-bulk expression profile : Hierarchical clustering of pseudo-bulk gene expression profiles (Wohnhaas 2019)


pilot.bulk_RNA_counts<-AverageExpression(pilot.integrated, group.by="Condition", slot= "counts",return.seurat = FALSE)
pilot.bulk_RNA_counts<-as.matrix(pilot.bulk_RNA_counts$RNA)

pilot.bulk_RNA_counts_CPM<-log2(cpm(pilot.bulk_RNA_counts, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.5))

pheatmap(pilot.bulk_RNA_counts_CPM,show_rownames=F,cluster_cols = FALSE)



#To further compare : calculate the pairwise correlation for all the cells profiled + clustering + heat map showing both (Fig. 1 D). (Wang 2021 )

pilot_counts_matrix<-as.matrix(pilot.integrated@assays[["RNA"]]@counts)

Pilot_counts_cor<-cor(pilot_counts_matrix)
Condition_df<-data.frame(pilot.integrated$Condition)
colnames(Condition_df)<-"Condition"
my_colour<-list(Condition=c(Fresh_cells = "#F8766D", Fixed_cells = "#00BFC4"))
pheatmap(Pilot_counts_cor,show_rownames=F,show_colnames=F,annotation_row =Condition_df,annotation_col =Condition_df,
         annotation_colors=my_colour,annotation_names_row = FALSE, annotation_names_col = FALSE,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")









#Now with cell identities
DefaultAssay(pilot.integrated) <- "integrated"
pilot.integrated <- RunPCA(pilot.integrated, features = VariableFeatures(object = pilot.integrated))
#Cluster the cells (12 dims)
pilot.integrated <- FindNeighbors(pilot.integrated, dims = 1:12, k.param=10)
pilot.integrated<- FindClusters(pilot.integrated, resolution = 0.5, algorithm=1)

#After cluster annotation with Linnarsson's database
#Assigning cell type identity to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("Astrocytes_RGL", "nIPCs_Neuroblasts", "Immature_neurons", 
                     "Mature_neurons","Microglia","Oligodendrocytes", "OPCs", "Endothelial")
pilot.integrated@active.ident <- plyr::mapvalues(x = pilot.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
pilot.integrated$CellType <- Idents(pilot.integrated)
#UMAP
pilot.integrated <- RunUMAP(pilot.integrated, dims = 1:15, n.neighbors = 15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pilot.integrated, reduction = "umap", label=TRUE, pt.size=2.3)



Distance_matrix<-dist(pilot.integrated@graphs[["integrated_nn"]]@p)



pilot_counts_matrix<-as.matrix(pilot.integrated@assays[["RNA"]]@counts)

Pilot_counts_cor<-cor(pilot_counts_matrix)
Condition_and_cluster_df<-data.frame(pilot.integrated$Condition)
Clusters<-as.factor(pilot.integrated$CellType)
Condition_and_cluster_df$Clusters<-Clusters
colnames(Condition_and_cluster_df)<-c("Condition", "Clusters")
my_colour<-list(Condition=c(Fresh_cells = "#F8766D", Fixed_cells = "#00BFC4"))
pheatmap(Pilot_counts_cor,show_rownames=F,show_colnames=F,annotation_row =Condition_and_cluster_df,annotation_col =Condition_and_cluster_df,
         annotation_colors=my_colour,annotation_names_row = FALSE, annotation_names_col = FALSE,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")


#Now I will do it with pilot.integrated@assays[["integrated"]]@scale.data


pilot_scaled_matrix<-as.matrix(pilot.integrated@assays[["integrated"]]@scale.data)

Pilot_scaled_cor<-cor(pilot_scaled_matrix)
Condition_and_cluster_df<-data.frame(pilot.integrated$Condition)
Clusters<-as.factor(pilot.integrated$CellType)
Condition_and_cluster_df$Clusters<-Clusters
colnames(Condition_and_cluster_df)<-c("Condition", "Clusters")
my_colour<-list(Condition=c(Fresh_cells = "#F8766D", Fixed_cells = "#00BFC4"),
                Clusters=c(Astrocytes_RGL= "#F8766D", nIPCs_Neuroblasts= "#CD9600", Immature_neurons= "#7CAE00", 
                            Mature_neurons= "#00BE67",Microglia= "#00BFC4",Oligodendrocytes="#00A9FF", OPCs="#C77CFF", Endothelial="#FF61CC"))
pheatmap(Pilot_scaled_cor,show_rownames=F,show_colnames=F,annotation_row =Condition_and_cluster_df,annotation_col =Condition_and_cluster_df,
         annotation_colors=my_colour,annotation_names_row = FALSE, annotation_names_col = FALSE,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")














#Subset cell identities
Astrocytes_RGL<-subset(pilot.integrated, subset = CellType=="Astrocytes +  RGL")
nIPCs_Neuroblasts<-subset(pilot.integrated, subset = CellType=="nIPCs + Neuroblasts")
Immature_neurons<-subset(pilot.integrated, subset = CellType=="Immature neurons")
Mature_neurons<-subset(pilot.integrated, subset = CellType=="Mature neurons")
Microglia<-subset(pilot.integrated, subset = CellType=="Microglia")
Oligodendrocytes<-subset(pilot.integrated, subset = CellType=="Oligodendrocytes")
OPCs<-subset(pilot.integrated, subset = CellType=="OPCs")
Endothelial_cells<-subset(pilot.integrated, subset = CellType=="Endothelial")

#number of cells per cell identity
table(Idents(pilot.integrated), pilot.integrated$Condition)

Astrocytes_RGL_matrix<-as.matrix(Astrocytes_RGL@assays[["RNA"]]@counts)

Astrocytes_RGL_cor<-cor(Astrocytes_RGL_matrix)
Condition_df_Astrocytes_RGL<-data.frame(Astrocytes_RGL$Condition)
colnames(Condition_df_Astrocytes_RGL)<-"Condition"
my_colour<-list(Condition=c(Fresh_cells = "#F8766D", Fixed_cells = "#00BFC4"))
pheatmap(Astrocytes_RGL_cor,show_rownames=F,show_colnames=F,annotation_row =Condition_df_Astrocytes_RGL,annotation_col =Condition_df_Astrocytes_RGL,
         annotation_colors=my_colour,annotation_names_row = FALSE, annotation_names_col = FALSE,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

nIPCs_Neuroblasts_matrix<-as.matrix(nIPCs_Neuroblasts@assays[["RNA"]]@counts)

nIPCs_Neuroblasts_cor<-cor(nIPCs_Neuroblasts_matrix)
Condition_df_nIPCs_Neuroblasts<-data.frame(nIPCs_Neuroblasts$Condition)
colnames(Condition_df_nIPCs_Neuroblasts)<-"Condition"
my_colour<-list(Condition=c(Fresh_cells = "#F8766D", Fixed_cells = "#00BFC4"))
pheatmap(nIPCs_Neuroblasts_cor,show_rownames=F,show_colnames=F,annotation_row =Condition_df_nIPCs_Neuroblasts,annotation_col =Condition_df_nIPCs_Neuroblasts,
         annotation_colors=my_colour,annotation_names_row = FALSE, annotation_names_col = FALSE,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
