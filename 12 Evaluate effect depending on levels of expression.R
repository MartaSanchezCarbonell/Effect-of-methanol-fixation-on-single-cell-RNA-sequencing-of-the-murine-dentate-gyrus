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
pilot.bulk_RNA_counts<-AverageExpression(pilot.integrated, group.by="Condition", slot= "counts")$RNA

#Save 
saveRDS(pilot.bulk_RNA_counts,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/pilot.bulk_RNA_counts.rds")
#Load 
pilot.bulk_RNA_counts<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/pilot.bulk_RNA_counts.rds")

pilot.bulk_RNA_counts<-as.data.frame(pilot.bulk_RNA_counts)


#--------Gene classification in 20 Groups separated by level of expression in Fresh cells. Each groups have the equal number of genes, except the last one that has 6 genes less-----------------------------------------------------

#I need to order the rows in ascending levels of expression in fresh cells
pilot.bulk_RNA_counts<-pilot.bulk_RNA_counts[order(pilot.bulk_RNA_counts$Fresh_cells),]

#Create a vector with repeated values grouping the genes in 20 groups (I have each gene 2 times because I have the condition Fixed or Fresh as a column)
Groups<-rep(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),each=774)
Groups<-head(Groups,-6)
pilot.bulk_RNA_counts$Groups<-Groups

#Create a vector with repeated values grouping the genes in 20 groups (I put the ranges as names now)
Range<-split(pilot.bulk_RNA_counts, pilot.bulk_RNA_counts$Groups)
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

Sorted_Head_and_tails$Fresh_cells

Groups_range<-rep(c("0-0","0-0.016","0.016-0.016","0.016-0.031","0.031-0.047","0.047-0.063","0.063-0.078","0.078-0.109","0.109-0.141","0.141-0.172","0.172-0.203","0.203-0.250","0.250-0.297","0.297-0.359","0.359-0.438","0.438-0.531","0.531-0.672","0.672-0.922","0.922-1.485","1.485-506.040"),each=774)
Groups_range<-head(Groups_range,-6)
pilot.bulk_RNA_counts$Groups_range<-Groups_range

#Save 
saveRDS(pilot.bulk_RNA_counts,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/pilot.bulk_RNA_counts_df.rds")
#Load 
pilot.bulk_RNA_counts<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/pilot.bulk_RNA_counts_df.rds")


#Separate the data frames in Fresh and Fixed
colnames(pilot.bulk_RNA_counts)
pilot.bulk_RNA_counts_Fresh<-pilot.bulk_RNA_counts[c("Fresh_cells","Groups","Groups_range")]
pilot.bulk_RNA_counts_Fixed<-pilot.bulk_RNA_counts[c("Fixed_cells","Groups","Groups_range")]

#Adding a column called condition with Fresh or Fixed and changeing col name form Fixed or fresh into Levels of expression
pilot.bulk_RNA_counts_Fresh$Condition<-rep("Fresh", times=15474)
colnames(pilot.bulk_RNA_counts_Fresh)[1] <- "Levels_of_expression"
pilot.bulk_RNA_counts_Fixed$Condition<-rep("Fixed", times=15474)
colnames(pilot.bulk_RNA_counts_Fixed)[1] <- "Levels_of_expression"

Effect_methanol_depending_expression<-rbind(pilot.bulk_RNA_counts_Fresh,pilot.bulk_RNA_counts_Fixed)

#Save
saveRDS(Effect_methanol_depending_expression, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Effect_methanol_depending_length_to_plot.rds")
#Load
Effect_methanol_depending_expression<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Effect_methanol_depending_length_to_plot.rds")

#Now plot the box plots
ggplot(Effect_methanol_depending_expression, aes(x=factor(Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Levels_of_expression_without outliers")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+ theme_classic()+
  ylim(0,5)

#Zoom-in in the low expression part of the graph
Zoom_in<-subset(Effect_methanol_depending_expression, Groups=="1" | Groups=="2" | Groups=="3" | Groups=="4" | Groups=="5")
ggplot(Zoom_in, aes(x=factor(Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Zoom_in_Levels_of_expression_without outliers")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+ theme_classic()+
  ylim(0,0.1)

#------------------------------Statistics-------------------------------------------
#statistical test to compare the average expression of the different gene groups(comparisons dos a dos, maybe bonferroni correction)
#Normality check with Q-Q plot (or quantile-quantile plot) 
#draws the correlation between a given sample and the normal distribution. A 45-degree reference line is also plotted.
ggqqplot(pilot.bulk_RNA_counts_Fresh$Levels_of_expression) + ggtitle("Fresh")
ggqqplot(pilot.bulk_RNA_counts_Fixed$Levels_of_expression) + ggtitle("Fixed")
#I will do Kruskal-Wallis Test in R
#Compute summary statistics by groups:
#Fresh
Stats_sum_Fresh<-group_by(pilot.bulk_RNA_counts_Fresh, Groups) %>%
  summarise(
    count = n(),
    mean = mean(Levels_of_expression, na.rm = TRUE),
    sd = sd(Levels_of_expression, na.rm = TRUE),
    median = median(Levels_of_expression, na.rm = TRUE),
    IQR = IQR(Levels_of_expression, na.rm = TRUE)
  )
Groups_range<-c("0-0","0-0.016","0.016-0.016","0.016-0.031","0.031-0.047","0.047-0.063","0.063-0.078","0.078-0.109","0.109-0.141","0.141-0.172","0.172-0.203","0.203-0.250","0.250-0.297","0.297-0.359","0.359-0.438","0.438-0.531","0.531-0.672","0.672-0.922","0.922-1.485","1.485-506.040")
Stats_sum_Fresh<-Stats_sum_Fresh[order(Stats_sum_Fresh$Groups),]
Stats_sum_Fresh$Groups_range<-Groups_range
# Export the summary table
time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Stats_sum_Fresh, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/",paste(time.name),"Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Fresh_statisics_sum", col.names = TRUE, row.names = TRUE, append = FALSE)

#Fixed
Stats_sum_Fixed<-group_by(pilot.bulk_RNA_counts_Fixed, Groups) %>%
  summarise(
    count = n(),
    mean = mean(Levels_of_expression, na.rm = TRUE),
    sd = sd(Levels_of_expression, na.rm = TRUE),
    median = median(Levels_of_expression, na.rm = TRUE),
    IQR = IQR(Levels_of_expression, na.rm = TRUE)
  )
Groups_range<-c("0-0","0-0.016","0.016-0.016","0.016-0.031","0.031-0.047","0.047-0.063","0.063-0.078","0.078-0.109","0.109-0.141","0.141-0.172","0.172-0.203","0.203-0.250","0.250-0.297","0.297-0.359","0.359-0.438","0.438-0.531","0.531-0.672","0.672-0.922","0.922-1.485","1.485-506.040")
Stats_sum_Fixed<-Stats_sum_Fixed[order(Stats_sum_Fixed$Groups),]
Stats_sum_Fixed$Groups_range<-Groups_range
# Export the summary table
write.xlsx(Stats_sum_Fixed, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/",paste(time.name),"Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Fixed_statisics_sum", col.names = TRUE, row.names = TRUE, append = T)


#Compute Kruskal-Wallis test
K<-kruskal.test(Levels_of_expression ~ interaction(Condition, Groups), data = Effect_methanol_depending_expression)
K1<-K$p.value

#Multiple pairwise-comparison between condition
#Statistic test comparing Fixed and Fresh in each length group
#To keep the script symple I change the name from the dataframe from "Sorted_Effect_methanol_depending_length_to_plot" to "df"
df<-Effect_methanol_depending_expression
p_values<-list()
Statistic_W<-list()
#1
W_test_t_1<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "1", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "1", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_1<-W_test_t_1$p.value
Statistic_W$W_1<-W_test_t_1$statistic
#2
W_test_t_2<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "2", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "2", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_2<-W_test_t_2$p.value
Statistic_W$W_2<-W_test_t_2$statistic
#3
W_test_t_3<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "3", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "3", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_3<-W_test_t_3$p.value
Statistic_W$W_3<-W_test_t_3$statistic
#4
W_test_t_4<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "4", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "4", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_4<-W_test_t_4$p.value
Statistic_W$W_4<-W_test_t_4$statistic
#5
W_test_t_5<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "5", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "5", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_5<-W_test_t_5$p.value
Statistic_W$W_5<-W_test_t_5$statistic
#6
W_test_t_6<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "6", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "6", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_6<-W_test_t_6$p.value
Statistic_W$W_6<-W_test_t_6$statistic
#7
W_test_t_7<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "7", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "7", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_7<-W_test_t_7$p.value
Statistic_W$W_7<-W_test_t_7$statistic
#8
W_test_t_8<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "8", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "8", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_8<-W_test_t_8$p.value
Statistic_W$W_8<-W_test_t_8$statistic
#9
W_test_t_9<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "9", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "9", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_9<-W_test_t_9$p.value
Statistic_W$W_9<-W_test_t_9$statistic
#10
W_test_t_10<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "10", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "10", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_10<-W_test_t_10$p.value
Statistic_W$W_10<-W_test_t_10$statistic
#11
W_test_t_11<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "11", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "11", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_11<-W_test_t_11$p.value
Statistic_W$W_11<-W_test_t_11$statistic
#12
W_test_t_12<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "12", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "12", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_12<-W_test_t_12$p.value
Statistic_W$W_12<-W_test_t_12$statistic
#13
W_test_t_13<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "13", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "13", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_13<-W_test_t_13$p.value
Statistic_W$W_13<-W_test_t_13$statistic
#14
W_test_t_14<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "14", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "14", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_14<-W_test_t_14$p.value
Statistic_W$W_14<-W_test_t_14$statistic
#15
W_test_t_15<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "15", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "15", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_15<-W_test_t_15$p.value
Statistic_W$W_15<-W_test_t_15$statistic
#16
W_test_t_16<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "16", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "16", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_16<-W_test_t_16$p.value
Statistic_W$W_16<-W_test_t_16$statistic
#17
W_test_t_17<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "17", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "17", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_17<-W_test_t_17$p.value
Statistic_W$W_17<-W_test_t_17$statistic
#18
W_test_t_18<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "18", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "18", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_18<-W_test_t_18$p.value
Statistic_W$W_18<-W_test_t_18$statistic
#19
W_test_t_19<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "19", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "19", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_19<-W_test_t_19$p.value
Statistic_W$W_19<-W_test_t_19$statistic
#20
W_test_t_20<-wilcox.test(df[df$Condition == "Fresh" & df$Groups == "20", ]$Levels_of_expression, df[df$Condition == "Fixed" & df$Groups == "20", ]$Levels_of_expression, paired=FALSE)
p_values$p_value_20<-W_test_t_20$p.value
Statistic_W$W_20<-W_test_t_20$statistic

#p-values adjustment
B_corrected_p_values<-as.list(p.adjust(p_values,method="bonferroni", n=length(p_values)))
B_Significative<-list()
for (i in 1:20){B_Significative[[i]]<-
  if (B_corrected_p_values[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values<-as.list(p.adjust(p_values,method="BH"), n=length(p_values))
BH_Significative<-list()
for (i in 1:20){BH_Significative[[i]]<-
  if (BH_corrected_p_values[[i]]<=0.05){"YES"} else {"NO"}}

p_values<-paste(p_values, collapse = NULL)
Statistic_W<-paste(Statistic_W, collapse = NULL)
B_corrected_p_values<-paste(B_corrected_p_values, collapse = NULL)
B_Significative<-paste(B_Significative, collapse = NULL)
BH_corrected_p_values<-paste(BH_corrected_p_values, collapse = NULL)
BH_Significative<-paste(BH_Significative, collapse = NULL)

#Results
Expression_groups<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")
Groups_range<-c("0-0","0-0.016","0.016-0.016","0.016-0.031","0.031-0.047","0.047-0.063","0.063-0.078","0.078-0.109","0.109-0.141","0.141-0.172","0.172-0.203","0.203-0.250","0.250-0.297","0.297-0.359","0.359-0.438","0.438-0.531","0.531-0.672","0.672-0.922","0.922-1.485","1.485-506.040")
Expression_groups_statistics_summary<-data.frame(Expression_groups,Groups_range,p_values,Statistic_W,B_corrected_p_values,B_Significative,BH_corrected_p_values,BH_Significative)
Expression_groups_statistics_summary

# Export the summary table
write.xlsx(Expression_groups_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/",paste(time.name),"Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Fresh_vs_Fixed", col.names = TRUE, row.names = TRUE, append = T)


#---------To really really see how is the distribution between conditions I need to:------------------------------------------------------------------------------
#Ratio calculation to see how much fixed cells decrease (ratio per gene group) + 
List_ratio<-split(Effect_methanol_depending_expression, Effect_methanol_depending_expression$Condition)

Ratio<-((1+(List_ratio[["Fresh"]][["Levels_of_expression"]]))/(1+(List_ratio[["Fixed"]][["Levels_of_expression"]])))
Ratio_df<-as.data.frame(Ratio)
Ratio_df$Groups<-List_ratio[["Fixed"]][["Groups"]]
Ratio_df$Groups_range<-List_ratio[["Fixed"]][["Groups_range"]]
#Save 
saveRDS(Ratio_df,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Ratio_df.rds")
#Load 
Ratio_df<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Ratio_df.rds")

ggplot(Ratio_df, aes(x=factor(Groups), y=Ratio)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Ratio (1+Fresh/1+Fixed) Average expression level as a function of transcript length")+
  scale_fill_manual(values=c("#9F9F9F"))+ 
  ylim(0.625,2)+ geom_hline(yintercept=1)+
  xlab("Expression categories") + ylab("Ratio 1+average expression Fresh/Fixed")+ theme_classic()


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
Groups_range<-c("0-0","0-0.016","0.016-0.016","0.016-0.031","0.031-0.047","0.047-0.063","0.063-0.078","0.078-0.109","0.109-0.141","0.141-0.172","0.172-0.203","0.203-0.250","0.250-0.297","0.297-0.359","0.359-0.438","0.438-0.531","0.531-0.672","0.672-0.922","0.922-1.485","1.485-506.040")
Stats_sum_Ratio<-Stats_sum_Ratio[order(Stats_sum_Ratio$Groups),]
Stats_sum_Ratio$Groups_range<-Groups_range
# Export the summary table
write.xlsx(Stats_sum_Ratio, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/",paste(time.name),"Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Ratio_statisics_sum", col.names = TRUE, row.names = TRUE, append = TRUE)

#Compute Kruskal-Wallis test
kruskal.test(Ratio ~ Groups, data = Ratio_df)

#Multiple pairwise-comparison between groups
#------Bonferroni correction
pairwise.wilcox.test_Bonferroni<-pairwise.wilcox.test(Ratio_df$Ratio, Ratio_df$Groups,
                                                      p.adjust.method = "bonferroni")
write.xlsx(pairwise.wilcox.test_Bonferroni$p.value, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/",paste(time.name),"Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Ratio_Bonferroni", col.names = TRUE, row.names = TRUE, append = TRUE)

#------FDR correction
pairwise.wilcox.test_BH<-pairwise.wilcox.test(Ratio_df$Ratio, Ratio_df$Groups,
                                              p.adjust.method = "BH")
write.xlsx(pairwise.wilcox.test_BH$p.value, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/",paste(time.name),"Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Ratio_Benjamini_Hochberg ", col.names = TRUE, row.names = TRUE, append = TRUE)


#To see if we obtain the same distribution with decontaminated counts
#I will add the two columns of decontaminated counts to the dataframe and check the expression levels, but keeping the original groups
#Therefore,
#Load
#Load 
pilot.bulk_RNA_counts<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/pilot.bulk_RNA_counts_df.rds")
pilot.bulk_RNA_counts$Genes<-row.names(pilot.bulk_RNA_counts)
#Load 
Corrected_Fresh<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/R objects/Corrected_Fresh.rds")
Corrected_Fixed<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/8 Ambient RNA investigation/R objects/Corrected_Fixed.rds")
#Create dataframes of the decontaminated counts
#RNA@counts
DefaultAssay(Corrected_Fresh) <- "RNA"
pilot.corrected_bulk_RNA_counts_Fresh<-as.data.frame(AverageExpression(Corrected_Fresh, slot= "counts")$RNA)
colnames(pilot.corrected_bulk_RNA_counts_Fresh)<-c("Decont.Fresh")
pilot.corrected_bulk_RNA_counts_Fresh$Genes<-row.names(pilot.corrected_bulk_RNA_counts_Fresh)

DefaultAssay(Corrected_Fixed) <- "RNA"
pilot.corrected_bulk_RNA_counts_Fixed<-as.data.frame(AverageExpression(Corrected_Fixed, slot= "counts")$RNA)
colnames(pilot.corrected_bulk_RNA_counts_Fixed)<-c("Decont.Fixed")
pilot.corrected_bulk_RNA_counts_Fixed$Genes<-row.names(pilot.corrected_bulk_RNA_counts_Fixed)

#Sort the dataframes by gene name
pilot.corrected_bulk_RNA_counts_Fresh<-pilot.corrected_bulk_RNA_counts_Fresh[order(pilot.corrected_bulk_RNA_counts_Fresh$Genes),]
pilot.corrected_bulk_RNA_counts_Fixed<-pilot.corrected_bulk_RNA_counts_Fixed[order(pilot.corrected_bulk_RNA_counts_Fixed$Genes),]
pilot.bulk_RNA_counts<-pilot.bulk_RNA_counts[order(pilot.bulk_RNA_counts$Genes),]

#Add decontaminated counts to "pilot.bulk_RNA_counts"
pilot.bulk_RNA_counts$Decont.Fresh<-pilot.corrected_bulk_RNA_counts_Fresh$Decont.Fresh
pilot.bulk_RNA_counts$Decont.Fixed<-pilot.corrected_bulk_RNA_counts_Fixed$Decont.Fixed

#Separate the data frames in Fresh and Fixed decontaminated counts
colnames(pilot.bulk_RNA_counts)
pilot.bulk_RNA_counts_Fresh<-pilot.bulk_RNA_counts[c("Decont.Fresh","Groups","Groups_range")]
pilot.bulk_RNA_counts_Fixed<-pilot.bulk_RNA_counts[c("Decont.Fixed","Groups","Groups_range")]

#Adding a column called condition with Fresh or Fixed and changeing col name form Fixed or fresh into Levels of expression
pilot.bulk_RNA_counts_Fresh$Condition<-rep("Fresh", times=15474)
colnames(pilot.bulk_RNA_counts_Fresh)[1] <- "Levels_of_expression"
pilot.bulk_RNA_counts_Fixed$Condition<-rep("Fixed", times=15474)
colnames(pilot.bulk_RNA_counts_Fixed)[1] <- "Levels_of_expression"

Effect_methanol_depending_expression_decontaminated<-rbind(pilot.bulk_RNA_counts_Fresh,pilot.bulk_RNA_counts_Fixed)

#Save
saveRDS(Effect_methanol_depending_expression_decontaminated, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Effect_methanol_depending_decontaminated_expression_to_plot.rds")
#Load
Effect_methanol_depending_expression_decontaminated<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Effect_methanol_depending_decontaminated_expression_to_plot.rds")

#Now plot the box plots
ggplot(Effect_methanol_depending_expression_decontaminated, aes(x=factor(Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Levels_of_decontaminated_expression_without outliers")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+ theme_classic()+
  ylim(0,5)

#Zoom-in in the low expression part of the graph
Zoom_in<-subset(Effect_methanol_depending_expression_decontaminated, Groups=="1" | Groups=="2" | Groups=="3" | Groups=="4" | Groups=="5")
ggplot(Zoom_in, aes(x=factor(Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Zoom_in_Levels_of_decontaminated_expression_without outliers")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+ theme_classic()+
  ylim(0,0.1)


#------------------------------Statistics-------------------------------------------
#statistical test to compare the average expression of the different gene groups(comparisons dos a dos, maybe bonferroni correction)
#Normality check with Q-Q plot (or quantile-quantile plot) 
#draws the correlation between a given sample and the normal distribution. A 45-degree reference line is also plotted.
ggqqplot(pilot.bulk_RNA_counts_Fresh$Levels_of_expression) + ggtitle("Fresh")
ggqqplot(pilot.bulk_RNA_counts_Fixed$Levels_of_expression) + ggtitle("Fixed")
#I will do Kruskal-Wallis Test in R
#Compute summary statistics by groups:
#Fresh
Stats_sum_Fresh<-group_by(pilot.bulk_RNA_counts_Fresh, Groups) %>%
  summarise(
    count = n(),
    mean = mean(Levels_of_expression, na.rm = TRUE),
    sd = sd(Levels_of_expression, na.rm = TRUE),
    median = median(Levels_of_expression, na.rm = TRUE),
    IQR = IQR(Levels_of_expression, na.rm = TRUE)
  )
Groups_range<-c("0-0","0-0.016","0.016-0.016","0.016-0.031","0.031-0.047","0.047-0.063","0.063-0.078","0.078-0.109","0.109-0.141","0.141-0.172","0.172-0.203","0.203-0.250","0.250-0.297","0.297-0.359","0.359-0.438","0.438-0.531","0.531-0.672","0.672-0.922","0.922-1.485","1.485-506.040")
Stats_sum_Fresh<-Stats_sum_Fresh[order(Stats_sum_Fresh$Groups),]
Stats_sum_Fresh$Groups_range<-Groups_range
# Export the summary table
write.xlsx(Stats_sum_Fresh, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Fresh_decont_statisics_sum", col.names = TRUE, row.names = TRUE, append = TRUE)
#Fixed
Stats_sum_Fixed<-group_by(pilot.bulk_RNA_counts_Fixed, Groups) %>%
  summarise(
    count = n(),
    mean = mean(Levels_of_expression, na.rm = TRUE),
    sd = sd(Levels_of_expression, na.rm = TRUE),
    median = median(Levels_of_expression, na.rm = TRUE),
    IQR = IQR(Levels_of_expression, na.rm = TRUE)
  )
Groups_range<-c("0-0","0-0.016","0.016-0.016","0.016-0.031","0.031-0.047","0.047-0.063","0.063-0.078","0.078-0.109","0.109-0.141","0.141-0.172","0.172-0.203","0.203-0.250","0.250-0.297","0.297-0.359","0.359-0.438","0.438-0.531","0.531-0.672","0.672-0.922","0.922-1.485","1.485-506.040")
Stats_sum_Fixed<-Stats_sum_Fixed[order(Stats_sum_Fixed$Groups),]
Stats_sum_Fixed$Groups_range<-Groups_range
# Export the summary table
write.xlsx(Stats_sum_Fixed, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Fixed_decont_statisics_sum", col.names = TRUE, row.names = TRUE, append = TRUE)

#---------To really really see how is the distribution between conditions I need to:------------------------------------------------------------------------------
#Ratio calculation to see how much fixed cells decrease (ratio per gene group) + 
List_Ratio_dec<-split(Effect_methanol_depending_expression_decontaminated, Effect_methanol_depending_expression_decontaminated$Condition)

Ratio_dec<-((1+(List_Ratio_dec[["Fresh"]][["Levels_of_expression"]]))/(1+(List_Ratio_dec[["Fixed"]][["Levels_of_expression"]])))
Ratio_dec_df<-as.data.frame(Ratio_dec)
Ratio_dec_df$Groups<-List_Ratio_dec[["Fixed"]][["Groups"]]
Ratio_dec_df$Groups_range<-List_Ratio_dec[["Fixed"]][["Groups_range"]]
#Save 
saveRDS(Ratio_dec_df,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Ratio_dec_df.rds")
#Load 
Ratio_dec_df<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Ratio_dec_df.rds")


ggplot(Ratio_dec_df, aes(x=factor(Groups), y=Ratio_dec)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Ratio_dec (1+Fresh/1+Fixed) Average expression level as a function of transcript length")+
  scale_fill_manual(values=c("#9F9F9F"))+ 
  ylim(0.625,2)+ geom_hline(yintercept=1)+
  xlab("Expression categories") + ylab("Ratio_dec 1+average expression Fresh/Fixed")+ theme_classic()


#statistical test to compare the Ratio_decs of the different gene groups(comparisons dos a dos, maybe bonferroni correction)
#Normality check with Q-Q plot (or quantile-quantile plot) 
#draws the correlation between a given sample and the normal distribution. A 45-degree reference line is also plotted.
ggqqplot(Ratio_dec_df$Ratio_dec)
#I will do Kruskal-Wallis Test in R
#Compute summary statistics by groups:
Stats_sum_Ratio_dec<-group_by(Ratio_dec_df, Groups) %>%
  summarise(
    count = n(),
    mean = mean(Ratio_dec, na.rm = TRUE),
    sd = sd(Ratio_dec, na.rm = TRUE),
    median = median(Ratio_dec, na.rm = TRUE),
    IQR = IQR(Ratio_dec, na.rm = TRUE)
  )
Groups_range<-c("0-0","0-0.016","0.016-0.016","0.016-0.031","0.031-0.047","0.047-0.063","0.063-0.078","0.078-0.109","0.109-0.141","0.141-0.172","0.172-0.203","0.203-0.250","0.250-0.297","0.297-0.359","0.359-0.438","0.438-0.531","0.531-0.672","0.672-0.922","0.922-1.485","1.485-506.040")
Stats_sum_Ratio_dec<-Stats_sum_Ratio_dec[order(Stats_sum_Ratio_dec$Groups),]
Stats_sum_Ratio_dec$Groups_range<-Groups_range
write.xlsx(Stats_sum_Ratio_dec, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Stats_sum_Ratio_decont", col.names = TRUE, row.names = TRUE, append = TRUE)


#-----------------------Multiple pairwise-comparison between ratio and decontaminated ratio----------------------------
#Statistic test comparing ratio and decontaminated ratio in each length group
#Load 
Ratio_df<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Ratio_df.rds")
#Load 
Ratio_dec_df<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/R objects/Ratio_dec_df.rds")


Ratio_df$Condition<-rep("Ratio", times=15474)
colnames(Ratio_df)[1] <- "Ratio"
Ratio_dec_df$Condition<-rep("Ratio_dec", times=15474)
colnames(Ratio_dec_df)[1] <- "Ratio"

Ratio_for_stat<-rbind(Ratio_df,Ratio_dec_df)

#To keep the script symple I change the name from the dataframe from "Sorted_Effect_methanol_depending_length_to_plot" to "df"
df<-Ratio_for_stat
p_values<-list()
Statistic_W<-list()
#1
W_test_t_1<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "1", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "1", ]$Ratio, paired=FALSE)
p_values$p_value_1<-W_test_t_1$p.value
Statistic_W$W_1<-W_test_t_1$statistic
#2
W_test_t_2<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "2", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "2", ]$Ratio, paired=FALSE)
p_values$p_value_2<-W_test_t_2$p.value
Statistic_W$W_2<-W_test_t_2$statistic
#3
W_test_t_3<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "3", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "3", ]$Ratio, paired=FALSE)
p_values$p_value_3<-W_test_t_3$p.value
Statistic_W$W_3<-W_test_t_3$statistic
#4
W_test_t_4<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "4", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "4", ]$Ratio, paired=FALSE)
p_values$p_value_4<-W_test_t_4$p.value
Statistic_W$W_4<-W_test_t_4$statistic
#5
W_test_t_5<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "5", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "5", ]$Ratio, paired=FALSE)
p_values$p_value_5<-W_test_t_5$p.value
Statistic_W$W_5<-W_test_t_5$statistic
#6
W_test_t_6<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "6", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "6", ]$Ratio, paired=FALSE)
p_values$p_value_6<-W_test_t_6$p.value
Statistic_W$W_6<-W_test_t_6$statistic
#7
W_test_t_7<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "7", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "7", ]$Ratio, paired=FALSE)
p_values$p_value_7<-W_test_t_7$p.value
Statistic_W$W_7<-W_test_t_7$statistic
#8
W_test_t_8<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "8", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "8", ]$Ratio, paired=FALSE)
p_values$p_value_8<-W_test_t_8$p.value
Statistic_W$W_8<-W_test_t_8$statistic
#9
W_test_t_9<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "9", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "9", ]$Ratio, paired=FALSE)
p_values$p_value_9<-W_test_t_9$p.value
Statistic_W$W_9<-W_test_t_9$statistic
#10
W_test_t_10<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "10", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "10", ]$Ratio, paired=FALSE)
p_values$p_value_10<-W_test_t_10$p.value
Statistic_W$W_10<-W_test_t_10$statistic
#11
W_test_t_11<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "11", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "11", ]$Ratio, paired=FALSE)
p_values$p_value_11<-W_test_t_11$p.value
Statistic_W$W_11<-W_test_t_11$statistic
#12
W_test_t_12<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "12", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "12", ]$Ratio, paired=FALSE)
p_values$p_value_12<-W_test_t_12$p.value
Statistic_W$W_12<-W_test_t_12$statistic
#13
W_test_t_13<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "13", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "13", ]$Ratio, paired=FALSE)
p_values$p_value_13<-W_test_t_13$p.value
Statistic_W$W_13<-W_test_t_13$statistic
#14
W_test_t_14<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "14", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "14", ]$Ratio, paired=FALSE)
p_values$p_value_14<-W_test_t_14$p.value
Statistic_W$W_14<-W_test_t_14$statistic
#15
W_test_t_15<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "15", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "15", ]$Ratio, paired=FALSE)
p_values$p_value_15<-W_test_t_15$p.value
Statistic_W$W_15<-W_test_t_15$statistic
#16
W_test_t_16<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "16", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "16", ]$Ratio, paired=FALSE)
p_values$p_value_16<-W_test_t_16$p.value
Statistic_W$W_16<-W_test_t_16$statistic
#17
W_test_t_17<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "17", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "17", ]$Ratio, paired=FALSE)
p_values$p_value_17<-W_test_t_17$p.value
Statistic_W$W_17<-W_test_t_17$statistic
#18
W_test_t_18<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "18", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "18", ]$Ratio, paired=FALSE)
p_values$p_value_18<-W_test_t_18$p.value
Statistic_W$W_18<-W_test_t_18$statistic
#19
W_test_t_19<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "19", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "19", ]$Ratio, paired=FALSE)
p_values$p_value_19<-W_test_t_19$p.value
Statistic_W$W_19<-W_test_t_19$statistic
#20
W_test_t_20<-wilcox.test(df[df$Condition == "Ratio" & df$Groups == "20", ]$Ratio, df[df$Condition == "Ratio_dec" & df$Groups == "20", ]$Ratio, paired=FALSE)
p_values$p_value_20<-W_test_t_20$p.value
Statistic_W$W_20<-W_test_t_20$statistic

#p-values adjustment
B_corrected_p_values<-as.list(p.adjust(p_values,method="bonferroni", n=length(p_values)))
B_Significative<-list()
for (i in 1:20){B_Significative[[i]]<-
  if (B_corrected_p_values[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values<-as.list(p.adjust(p_values,method="BH"), n=length(p_values))
BH_Significative<-list()
for (i in 1:20){BH_Significative[[i]]<-
  if (BH_corrected_p_values[[i]]<=0.05){"YES"} else {"NO"}}

p_values<-paste(p_values, collapse = NULL)
Statistic_W<-paste(Statistic_W, collapse = NULL)
B_corrected_p_values<-paste(B_corrected_p_values, collapse = NULL)
B_Significative<-paste(B_Significative, collapse = NULL)
BH_corrected_p_values<-paste(BH_corrected_p_values, collapse = NULL)
BH_Significative<-paste(BH_Significative, collapse = NULL)

#Results
Expression_groups<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")
Groups_range<-c("0-0","0-0.016","0.016-0.016","0.016-0.031","0.031-0.047","0.047-0.063","0.063-0.078","0.078-0.109","0.109-0.141","0.141-0.172","0.172-0.203","0.203-0.250","0.250-0.297","0.297-0.359","0.359-0.438","0.438-0.531","0.531-0.672","0.672-0.922","0.922-1.485","1.485-506.040")
Ratio_groups_statistics_summary<-data.frame(Expression_groups,Groups_range,p_values,Statistic_W,B_corrected_p_values,B_Significative,BH_corrected_p_values,BH_Significative)
Ratio_groups_statistics_summary

# Export the summary table
write.xlsx(Ratio_groups_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/12 Evaluate effect depending on levels of expression/Effect_expression_levels_statistics_summary.xlsx"))
           , sheetName = "Ratio_vs_Ratio_dec", col.names = TRUE, row.names = TRUE, append = T)
