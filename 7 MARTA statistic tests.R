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
install.packages("ggdendro")
devtools::install_github("YuLab-SMU/ggtree")
install.packages("patchwork")
install.packages("readxl")
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
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 
library(readxl)


#Load the data
Fresh_cells<-read_excel("F:\\Marta\\Getting started with scRNAseq analysis\\Pilot results\\Pilot analysis + results\\7 Evaluate dropouts PATRICIA\\Statistics Dropouts.xlsx"
                        ,sheet = "Fresh")
Fixed_cells<-read_excel("F:\\Marta\\Getting started with scRNAseq analysis\\Pilot results\\Pilot analysis + results\\7 Evaluate dropouts PATRICIA\\Statistics Dropouts.xlsx"
                        ,sheet = "Fixed")

W_test_t_0.05<-wilcox.test(Fresh_cells$`0.05`, Fixed_cells$`0.05`, paired=FALSE)
p_value_0.05<-W_test_t_0.05$p.value
W_0.05<-W_test_t_0.05$statistic
B_corrected_p_value_0.05<-p_value_0.05*5
Significative_0.05_B<-if(B_corrected_p_value_0.05<=0.05) {
  "YES"
} else {
  "NO"
}

W_test_t_2<-wilcox.test(Fresh_cells$`2`, Fixed_cells$`2`, paired=FALSE)
p_value_2<-W_test_t_2$p.value
W_2<-W_test_t_2$statistic
B_corrected_p_value_2<-p_value_2*5
Significative_2_B<-if(B_corrected_p_value_2<=0.05) {
  "YES"
} else {
  "NO"
}

W_test_t_6<-wilcox.test(Fresh_cells$`6`, Fixed_cells$`6`, paired=FALSE)
p_value_6<-W_test_t_6$p.value
W_6<-W_test_t_6$statistic
B_corrected_p_value_6<-p_value_6*5
Significative_6_B<-if(B_corrected_p_value_6<=0.05) {
  "YES"
} else {
  "NO"
}

W_test_t_10<-wilcox.test(Fresh_cells$`10`, Fixed_cells$`10`, paired=FALSE)
p_value_10<-W_test_t_10$p.value
W_10<-W_test_t_10$statistic
B_corrected_p_value_10<-p_value_10*5
Significative_10_B<-if(B_corrected_p_value_10<=0.05) {
  "YES"
} else {
  "NO"
}

W_test_t_25<-wilcox.test(Fresh_cells$`25`, Fixed_cells$`25`, paired=FALSE)
p_value_25<-W_test_t_25$p.value
W_25<-W_test_t_25$statistic
B_corrected_p_value_25<-p_value_25*5
Significative_25_B<-if(B_corrected_p_value_25<=0.05) {
  "YES"
} else {
  "NO"
}

#p-values adjustment
p_values<-c(p_value_0.05,p_value_2,p_value_6,p_value_10,p_value_25)
B_corrected_p_values<-p.adjust(p_values,method="bonferroni")
BH_corrected_p_values<-p.adjust(p_values,method="BH")


#Results
Thresholds<-c("0.05", "2", "6","10","25")
p_value<-c(p_value_0.05,p_value_2,p_value_6,p_value_10,p_value_25)
Statistic_W<-c(W_0.05,W_2,W_6,W_10,W_25)
Threshholds_statistics_summary<-data.frame(Thresholds,p_value,Statistic_W,B_corrected_p_values,BH_corrected_p_values)
Threshholds_statistics_summary
