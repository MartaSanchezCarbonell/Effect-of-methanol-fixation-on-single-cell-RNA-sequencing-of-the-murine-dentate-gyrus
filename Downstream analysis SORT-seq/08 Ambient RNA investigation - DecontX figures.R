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
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
BiocManager::install("celda")
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
library(SingleCellExperiment)
library(forcats)
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


#Loading into Seurat
Corrected_Fresh= CreateSeuratObject(decontXcounts(sce.Fresh_cells))
Corrected_Fixed= CreateSeuratObject(decontXcounts(sce.Fixed_cells))
Corrected_Fresh[["Contamination_fraction"]]<-sce.Fresh_cells@metadata[["decontX"]][["contamination"]]
Corrected_Fixed[["Contamination_fraction"]]<-sce.Fixed_cells@metadata[["decontX"]][["contamination"]]


#Save 
saveRDS(Corrected_Fresh,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/R objects/Corrected_Fresh.rds")
saveRDS(Corrected_Fixed,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/R objects/Corrected_Fixed.rds")

#Load 
Corrected_Fresh<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/R objects/Corrected_Fresh.rds")
Corrected_Fixed<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/R objects/Corrected_Fixed.rds")




#Density plots with box plots
#Arrangement of the dataframe
Contamination_Fresh<-Corrected_Fresh$Contamination_fraction
Condition<-rep("Faresh", times=64)
Condition_num<-rep("1", times=64)
Contamination_Fresh_df<-data.frame(Contamination_Fresh,Condition,Condition_num)
colnames(Contamination_Fresh_df)[1] <- "Contamination"

Contamination_Fixed<-Corrected_Fixed$Contamination_fraction
Condition<-rep("Fixed", times=132)
Condition_num<-rep("2", times=132)
Contamination_Fixed_df<-data.frame(Contamination_Fixed,Condition,Condition_num)
colnames(Contamination_Fixed_df)[1] <- "Contamination"

Contamination_df<-rbind(Contamination_Fresh_df,Contamination_Fixed_df)
Contamination_df$Parameter<-rep("Contamination_fraction", times=196)

#Compute summary statistics by groups:
Stats_sum_Contamination<-group_by(Contamination_df, Condition) %>%
  summarise(
    count = n(),
    mean = mean(Contamination, na.rm = TRUE),
    sd = sd(Contamination, na.rm = TRUE),
    median = median(Contamination, na.rm = TRUE),
    IQR = IQR(Contamination, na.rm = TRUE)
  )

# Export the summary table
time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Stats_sum_Contamination, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/",paste(time.name),"Contamination_statistics_summary.xlsx"))
           , sheetName = "Contamination_statisics_sum", col.names = TRUE, row.names = TRUE, append = FALSE)

#Normality test
shapiro.test(Contamination_Fresh_df$Contamination)
shapiro.test(Contamination_Fixed_df$Contamination)

#Wilcoxon rank-Sum test

#Fresh Cells Before - Fixed cells Before
wilcox.test(Contamination_Fresh_df$Contamination, Contamination_Fixed_df$Contamination, paired=FALSE)


ggplot(Contamination_df, aes(x=Parameter, y=Contamination, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F) 







