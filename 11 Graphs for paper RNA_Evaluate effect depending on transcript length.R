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
BiocManager::install("ExperimentHub")
install.packages("pheatmap")
install.packages("Hmisc")
install.packages("Rtools")
#remove.packages("cli")
install.packages("cli")
BiocManager::install("GenomicFeatures")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Mm.eg.db")
install.packages('data.table')
install.packages("xlsx")
#Activate the following libraries/packages
library(rlang)
library(dplyr)
library(cli)
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
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(data.table)
library(xlsx)

#Load
Effect_methanol_depending_length<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Effect_methanol_depending_length.rds")


#I will try to plot all together
Effect_methanol_depending_length_Fresh<-Effect_methanol_depending_length[,-3]
Effect_methanol_depending_length_Fixed<-Effect_methanol_depending_length[,-2]
#Adding a column called condition with Fresh or Fixed and changeing col name form Fixed or fres into Levels of expression
Effect_methanol_depending_length_Fresh$Condition<-rep("Faresh", times=14100)
colnames(Effect_methanol_depending_length_Fresh)[2] <- "Levels_of_expression"
Effect_methanol_depending_length_Fixed$Condition<-rep("Fixed", times=14100)
colnames(Effect_methanol_depending_length_Fixed)[2] <- "Levels_of_expression"

Effect_methanol_depending_length_to_plot<-rbind(Effect_methanol_depending_length_Fresh,Effect_methanol_depending_length_Fixed)

#-------------------Gene classification in 20 Groups separated by length. Each groups have the equal number of genes-----------------------------------------------------

#I need to order the rows in ascending exonic.gene.sizes
Sorted_Effect_methanol_depending_length_to_plot<-Effect_methanol_depending_length_to_plot[order(Effect_methanol_depending_length_to_plot$exonic.gene.sizes),]

#Create a vector with repeated values grouping the genes in 20 groups (I have each gene 2 times because I have the condition Fixed or Fresh as a column)
Groups<-rep(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),each=1410)
Sorted_Effect_methanol_depending_length_to_plot$Groups<-Groups

#Save
saveRDS(Sorted_Effect_methanol_depending_length_to_plot, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Graphs_for_paper_Sorted_Effect_methanol_depending_length_to_plot.rds")
#Load
Sorted_Effect_methanol_depending_length_to_plot<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Graphs_for_paper_Sorted_Effect_methanol_depending_length_to_plot.rds")

#Box plots

ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=factor(Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Levels_of_expression_ without outliers")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"))+ 
  ylim(0,2)+ theme_classic()

#Plots to represent statistics
#Create a vector with repeated values grouping the genes in 20 groups (I put the ranges as names now)
Range<-split(Sorted_Effect_methanol_depending_length_to_plot, Sorted_Effect_methanol_depending_length_to_plot$Groups)
Heads<-data.frame()
for (i in 1:20){tmp <- head(Range[[i]], n=1)
Heads<- rbind(Heads, tmp)
}
Tails<-data.frame()
for (i in 1:20){tmp <- tail(Range[[i]], n=1)
Tails<- rbind(Tails, tmp)
}
Head_and_tails<-rbind(Heads,Tails)
Sorted_Head_and_tails<-Head_and_tails[order(Head_and_tails$Groups),]

Sorted_Head_and_tails$exonic.gene.sizes


Groups_range<-rep(c("66-1272","1272-1701","1702-2083","2084-2419","2419-2724","2724-3031","3033-3341","3342-3652","3652-3995","3995-4358","4358-4733","4734-5193","5194-5631","5633-6129","6130-6772","6774-7556","7556-8514","8514-9857","9859-12513","12520-117334"),each=1410)
Sorted_Effect_methanol_depending_length_to_plot$Groups_range<-Groups_range
#Comparison between Fresh and Fixed
#Save
saveRDS(Sorted_Effect_methanol_depending_length_to_plot, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Graphs_for_paper_Sorted_Effect_methanol_depending_length_to_plot_2.rds")
#Load
Sorted_Effect_methanol_depending_length_to_plot<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Graphs_for_paper_Sorted_Effect_methanol_depending_length_to_plot_2.rds")

ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=reorder(Groups_range,Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Average expression level as a function of transcript length")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"))+ 
  ylim(0,1.7)+
  xlab("Transcript length categories") + ylab("Average expression")+ theme_classic()


#-------------------------------------------------Statistics with Kruskal-Wallis Test in R for Fresh---------------------------------------------
#statistical test to compare the expression of the different gene groups(comparisons dos a dos, maybe bonferroni correction)
#Load
Sorted_Effect_methanol_depending_length_to_plot<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Sorted_Effect_methanol_depending_length_to_plot_2.rds")
Fresh_and_Fixed<-split(Sorted_Effect_methanol_depending_length_to_plot, Sorted_Effect_methanol_depending_length_to_plot$Condition)
Fresh<-Fresh_and_Fixed$Fresh
#Normality check with Q-Q plot (or quantile-quantile plot) 
#draws the correlation between a given sample and the normal distribution. A 45-degree reference line is also plotted.
ggqqplot(Fresh$Levels_of_expression)
#I will do Kruskal-Wallis Test in R
#Compute summary statistics by groups:
Stats_sum_Fresh<-group_by(Fresh, Groups) %>%
  summarise(
    count = n(),
    mean = mean(Levels_of_expression, na.rm = TRUE),
    sd = sd(Levels_of_expression, na.rm = TRUE),
    median = median(Levels_of_expression, na.rm = TRUE),
    IQR = IQR(Levels_of_expression, na.rm = TRUE)
  )
Group_Range<-c("66-1272","1272-1701","1702-2083","2084-2419","2419-2724","2724-3031","3033-3341","3342-3652","3652-3995","3995-4358","4358-4733","4734-5193","5194-5631","5633-6129","6130-6772","6774-7556","7556-8514","8514-9857","9859-12513","12520-117334")
Stats_sum_Fresh<-Stats_sum_Fresh[order(Stats_sum_Fresh$Groups),]
Stats_sum_Fresh$Groups_range<-Group_Range
# Export the summary table
time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Stats_sum_Fresh, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary_v02.xlsx"))
           , sheetName = "Fresh_statisics_sum", col.names = TRUE, row.names = TRUE, append = FALSE)

#Compute Kruskal-Wallis test
kruskal.test(Levels_of_expression ~ Groups, data = Fresh)

#Multiple pairwise-comparison between groups
#------Bonferroni correction
pairwise.wilcox.test_Bonferroni<-pairwise.wilcox.test(Fresh$Levels_of_expression, Fresh$Groups,
                                                      p.adjust.method = "bonferroni")
write.xlsx(pairwise.wilcox.test_Bonferroni$p.value, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary_v02.xlsx"))
           , sheetName = "Fresh_Bonferroni", col.names = TRUE, row.names = TRUE, append = TRUE)

#------FDR correction
pairwise.wilcox.test_BH<-pairwise.wilcox.test(Fresh$Levels_of_expression, Fresh$Groups,
                                              p.adjust.method = "BH")
write.xlsx(pairwise.wilcox.test_BH$p.value, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary_v02.xlsx"))
           , sheetName = "Fresh_Benjamini_Hochberg ", col.names = TRUE, row.names = TRUE, append = TRUE)

#-------------------------------------------------Statistics with Kruskal-Wallis Test in R for Fixed---------------------------------------------
#statistical test to compare the expression of the different gene groups(comparisons dos a dos, maybe bonferroni correction)
#Load
Sorted_Effect_methanol_depending_length_to_plot<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Sorted_Effect_methanol_depending_length_to_plot_2.rds")
Fresh_and_Fixed<-split(Sorted_Effect_methanol_depending_length_to_plot, Sorted_Effect_methanol_depending_length_to_plot$Condition)
Fixed<-Fresh_and_Fixed$Fixed
#Normality check with Q-Q plot (or quantile-quantile plot) 
#draws the correlation between a given sample and the normal distribution. A 45-degree reference line is also plotted.
ggqqplot(Fixed$Levels_of_expression)
#I will do Kruskal-Wallis Test in R
#Compute summary statistics by groups:
Stats_sum_Fixed<-group_by(Fixed, Groups) %>%
  summarise(
    count = n(),
    mean = mean(Levels_of_expression, na.rm = TRUE),
    sd = sd(Levels_of_expression, na.rm = TRUE),
    median = median(Levels_of_expression, na.rm = TRUE),
    IQR = IQR(Levels_of_expression, na.rm = TRUE)
  )
Group_Range<-c("66-1272","1272-1701","1702-2083","2084-2419","2419-2724","2724-3031","3033-3341","3342-3652","3652-3995","3995-4358","4358-4733","4734-5193","5194-5631","5633-6129","6130-6772","6774-7556","7556-8514","8514-9857","9859-12513","12520-117334")
Stats_sum_Fixed<-Stats_sum_Fixed[order(Stats_sum_Fixed$Groups),]
Stats_sum_Fixed$Groups_range<-Group_Range

write.xlsx(Stats_sum_Fixed, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary_v02.xlsx"))
           , sheetName = "Fixed_statisics_sum", col.names = TRUE, row.names = TRUE, append = T)

#Compute Kruskal-Wallis test
kruskal.test(Levels_of_expression ~ Groups, data = Fixed)

#Multiple pairwise-comparison between groups
#------Bonferroni correction
pairwise.wilcox.test_Bonferroni<-pairwise.wilcox.test(Fixed$Levels_of_expression, Fixed$Groups,
                                                      p.adjust.method = "bonferroni")
write.xlsx(pairwise.wilcox.test_Bonferroni$p.value, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary_v02.xlsx"))
           , sheetName = "Fixed_Bonferroni", col.names = TRUE, row.names = TRUE, append = TRUE)

#------FDR correction
pairwise.wilcox.test_BH<-pairwise.wilcox.test(Fixed$Levels_of_expression, Fixed$Groups,
                                              p.adjust.method = "BH")
write.xlsx(pairwise.wilcox.test_BH$p.value, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary_v02.xlsx"))
           , sheetName = "Fixed_Benjamini_Hochberg ", col.names = TRUE, row.names = TRUE, append = TRUE)



#--------------------------------------------------------To really really see if the distribution between conditions is similar I need to:------------------------------------------------------------------------------
#Ratio calculation to see how much fixed cells decrease (ratio per gene group) + 
List_ratio<-split(Sorted_Effect_methanol_depending_length_to_plot, Sorted_Effect_methanol_depending_length_to_plot$Condition)

Ratio<-((1+(List_ratio[["Fresh"]][["Levels_of_expression"]]))/(1+(List_ratio[["Fixed"]][["Levels_of_expression"]])))
Ratio_df<-as.data.frame(Ratio)
Ratio_df$SYMBOL<-List_ratio[["Fixed"]][["SYMBOL"]]
Ratio_df$ENSEMBL<-List_ratio[["Fixed"]][["ENSEMBL"]]
Ratio_df$exonic.gene.sizes<-List_ratio[["Fixed"]][["exonic.gene.sizes"]]
Ratio_df$Groups<-List_ratio[["Fixed"]][["Groups"]]
Ratio_df$Groups_range<-List_ratio[["Fixed"]][["Groups_range"]]

ggplot(Ratio_df, aes(x=factor(Groups), y=Ratio)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Ratio (1+Fresh/1+Fixed) Average expression level as a function of transcript length")+
  scale_fill_manual(values=c("#9F9F9F"))+ 
  ylim(0.625,1.75)+ geom_hline(yintercept=1)+
  xlab("Transcript length categories") + ylab("Ratio 1+average expression Fresh/Fixed") + theme_classic()


#statistical test to compare the ratios of the different gene groups(comparisons dos a dos, maybe bonferroni correction)
#Normality check with Q-Q plot (or quantile-quantile plot) 
#draws the correlation between a given sample and the normal distribution. A 45-degree reference line is also plotted.
ggqqplot(Ratio_df$Ratio)
#I will do Kruskal-Wallis Test in R
#Compute summary statistics by groups:
Stats_sum_Ratio<-group_by(Ratio_df, Groups) %>%
  summarise(
    count = n(),
    mean = mean(Ratio, na.rm = TRUE),
    sd = sd(Ratio, na.rm = TRUE),
    median = median(Ratio, na.rm = TRUE),
    IQR = IQR(Ratio, na.rm = TRUE)
  )
Group_Range<-c("66-1272","1272-1701","1702-2083","2084-2419","2419-2724","2724-3031","3033-3341","3342-3652","3652-3995","3995-4358","4358-4733","4734-5193","5194-5631","5633-6129","6130-6772","6774-7556","7556-8514","8514-9857","9859-12513","12520-117334")
Stats_sum_Ratio<-Stats_sum_Ratio[order(Stats_sum_Ratio$Groups),]
Stats_sum_Ratio$Groups_range<-Group_Range
# Export the summary table
time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Stats_sum_Ratio, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_RATIO_statistics_summary.xlsx"))
           , sheetName = "Ratio_statisics_sum", col.names = TRUE, row.names = TRUE, append = FALSE)

#Compute Kruskal-Wallis test
kruskal.test(Ratio ~ Groups, data = Ratio_df)

#Multiple pairwise-comparison between groups
#------Bonferroni correction
pairwise.wilcox.test_Bonferroni<-pairwise.wilcox.test(Ratio_df$Ratio, Ratio_df$Groups,
                     p.adjust.method = "bonferroni")
write.xlsx(pairwise.wilcox.test_Bonferroni$p.value, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_RATIO_statistics_summary.xlsx"))
           , sheetName = "Bonferroni", col.names = TRUE, row.names = TRUE, append = TRUE)

#------FDR correction
pairwise.wilcox.test_BH<-pairwise.wilcox.test(Ratio_df$Ratio, Ratio_df$Groups,
                     p.adjust.method = "BH")
write.xlsx(pairwise.wilcox.test_BH$p.value, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_RATIO_statistics_summary.xlsx"))
           , sheetName = "Benjamini_Hochberg ", col.names = TRUE, row.names = TRUE, append = TRUE)



