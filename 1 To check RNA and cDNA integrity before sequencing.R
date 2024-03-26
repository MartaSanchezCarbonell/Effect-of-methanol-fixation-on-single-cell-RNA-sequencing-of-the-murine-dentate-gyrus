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
library(celda)
library(vioplot)



Fresh_cells<-c(8.90,8.60,8.30,8.77,8.83)
Fixed_cells<-c(7.30,8.70,7.50,9.03,8.40)



#STATISTICS

#Means
Fresh_cells_mean<-mean(Fresh_cells)
Fresh_cells_mean
Fixed_cells_mean<-mean(Fixed_cells)
Fixed_cells_mean


#Median
Fresh_cells_median<-median(Fresh_cells)
Fresh_cells_median
Fixed_cells_median<-median(Fixed_cells)
Fixed_cells_median


#standard error of the mean (SEM)

std <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

Fresh_cells_std<-std(Fresh_cells)
Fresh_cells_std
Fixed_cells_std<-std(Fixed_cells)
Fixed_cells_std


#Summary
Condition<-c("Fresh cells","Fixed cells")
RIN_mean<-c(Fresh_cells_mean,Fixed_cells_mean)
RIN_median<-c(Fresh_cells_median,Fixed_cells_median)
RIN_std<-c(Fresh_cells_std,Fixed_cells_std)
Summary<-data.frame(Condition,RIN_mean,RIN_median,RIN_std)
Summary

#Normality test
shapiro.test(Fresh_cells)
shapiro.test(Fixed_cells)
#Homogeneity of Variance
var.test(Fresh_cells,Fixed_cells,alternative = "two.sided")

#Wilcox
wilcox.test(Fresh_cells,Fixed_cells, paired=TRUE)
#T-Test
t.test(Fresh_cells,Fixed_cells, paired=TRUE)