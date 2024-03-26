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
install.packages("viridis")
install.packages("hrbrthemes")

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
library(hrbrthemes)
library(viridis)


Methanol_Treh_recovery<-c(15.27,15.26,25.26,25.26,36.32,22.24,28.29,31.05,13.29,35.26,31.58,64.47)
Cryo_recovery<-c(5.56,2.61,2.13,6.17,2.69,3.81,4.13,21.30)
Methanol_recovery<-c(11.75, 11.45, 12.63, 11.32, 41.97, 77.50, 53.82)

#STATISTICS

#Means
Methanol_Treh_recovery_mean<-mean(Methanol_Treh_recovery)
Methanol_Treh_recovery_mean
Cryo_recovery_mean<-mean(Cryo_recovery)
Cryo_recovery_mean
Methanol_recovery_mean<-mean(Methanol_recovery)
Methanol_recovery_mean

#Median
Methanol_Treh_recovery_median<-median(Methanol_Treh_recovery)
Methanol_Treh_recovery_median
Cryo_recovery_median<-median(Cryo_recovery)
Cryo_recovery_median
Methanol_recovery_median<-median(Methanol_recovery)
Methanol_recovery_median

#standard error of the mean (SEM)

std <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

Methanol_Treh_recovery_std<-std(Methanol_Treh_recovery)
Methanol_Treh_recovery_std
Cryo_recovery_std<-std(Cryo_recovery)
Cryo_recovery_std
Methanol_recovery_std<-std(Methanol_recovery)
Methanol_recovery_std

#Summary
Preservation<-c("Methanol Treh","Cryo", "Methanol")
Recovery_rate_mean<-c(Methanol_Treh_recovery_mean,Cryo_recovery_mean, Methanol_recovery_mean)
Recovery_rate_median<-c(Methanol_Treh_recovery_median,Cryo_recovery_median, Methanol_recovery_median)
Recovery_rate_std<-c(Methanol_Treh_recovery_std,Cryo_recovery_std,Methanol_recovery_std )
Summary<-data.frame(Preservation,Recovery_rate_mean,Recovery_rate_median,Recovery_rate_std)
Summary


shapiro.test(Methanol_Treh_recovery)
shapiro.test(Cryo_recovery)
shapiro.test(Methanol_recovery)

wilcox.test(Methanol_Treh_recovery,Cryo_recovery, paired=FALSE)
wilcox.test(Methanol_recovery,Cryo_recovery, paired=FALSE)
wilcox.test(Methanol_Treh_recovery, Methanol_recovery, paired=FALSE)

#box plots
Preservations<-data.frame(
      method=c(rep("Methanol_Treh_recovery",12),rep("Methanol_Treh_recovery",12), rep("Cryo_recovery",8)),
      recovery=c(Methanol_Treh_recovery,Cryo_recovery)
)

ggplot(Preservations, aes(x=method, y=recovery,fill=method)) + 
  ylab("Recovery (%)")+
  geom_boxplot(outlier.shape=16)+ stat_summary(fun.y=mean, geom="point", shape=1, size=4)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_x_discrete(limits=c("Cryo_recovery","Methanol_Treh_recovery"))+
  scale_fill_manual(values=c("#FFFFFF","#00BFC4"))+
  theme_classic()+
  theme(legend.position="none",axis.title.x = element_blank())

#box plots without treh included
Preservations<-data.frame(
  method=c(rep("Methanol_recovery",7), rep("Methanol_Treh_recovery",12)),
  recovery=c(Methanol_recovery,Methanol_Treh_recovery)
)

ggplot(Preservations, aes(x=method, y=recovery,fill=method)) + 
  ylab("Recovery (%)")+
  geom_boxplot(outlier.shape=16)+ stat_summary(fun.y=mean, geom="point", shape=1, size=4)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_x_discrete(limits=c("Methanol_recovery","Methanol_Treh_recovery"))+
  scale_fill_manual(values=c("#FFFFFF","#00BFC4"))+
  theme_classic()+
  theme(legend.position="none",axis.title.x = element_blank())


