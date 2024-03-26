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
saveRDS(pilot.bulk_RNA_counts,file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/pilot.bulk_RNA_counts.rds")

#Load 
pilot.bulk_RNA_counts<-readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/pilot.bulk_RNA_counts.rds")



# First, import the GTF-file that you have also used as input for htseq-count
txdb <- makeTxDbFromGFF("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/Mus_musculus.GRCm38.102.gtf",format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
exonic.gene.sizes.dataframe<-as.data.frame(exonic.gene.sizes)
exonic.gene.sizes.dataframe$ENSEMBL<-row.names(exonic.gene.sizes.dataframe)
#Annotate the gene symbols
ENSEMBL_gene_names<-as.vector(row.names(exonic.gene.sizes.dataframe))
annots <- AnnotationDbi::select(org.Mm.eg.db, keys = ENSEMBL_gene_names,
                                columns = "SYMBOL", keytype ="ENSEMBL")
#merging annots with exonic.gene.sizes.dataframe
Total_Non_Overlapping_Exon_Length_Per_Gene<-merge(exonic.gene.sizes.dataframe, annots,by="ENSEMBL")

#Save
saveRDS(Total_Non_Overlapping_Exon_Length_Per_Gene, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Total_Non_Overlapping_Exon_Length_Per_Gene.rds")
#Load
Total_Non_Overlapping_Exon_Length_Per_Gene<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Total_Non_Overlapping_Exon_Length_Per_Gene.rds")

#Merging Total_Non_Overlapping_Exon_Length_Per_Gene with pilot.bulk_RNA_counts
pilot.bulk_RNA_counts_df<-as.data.frame(pilot.bulk_RNA_counts)
pilot.bulk_RNA_counts_df$SYMBOL<-row.names(pilot.bulk_RNA_counts_df)

Effect_methanol_depending_length<-merge(pilot.bulk_RNA_counts_df,Total_Non_Overlapping_Exon_Length_Per_Gene, by="SYMBOL")

#Save
saveRDS(Effect_methanol_depending_length, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Effect_methanol_depending_length.rds")
#Load
Effect_methanol_depending_length<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Effect_methanol_depending_length.rds")


#I will try to plot all together
Effect_methanol_depending_length_Fresh<-Effect_methanol_depending_length[,-3]
Effect_methanol_depending_length_Fixed<-Effect_methanol_depending_length[,-2]
#Adding a column called condition with Fresh or Fixed and changeing col name form Fixed or fres into Levels of expression
Effect_methanol_depending_length_Fresh$Condition<-rep("Fresh", times=14100)
colnames(Effect_methanol_depending_length_Fresh)[2] <- "Levels_of_expression"
Effect_methanol_depending_length_Fixed$Condition<-rep("Fixed", times=14100)
colnames(Effect_methanol_depending_length_Fixed)[2] <- "Levels_of_expression"

Effect_methanol_depending_length_to_plot<-rbind(Effect_methanol_depending_length_Fresh,Effect_methanol_depending_length_Fixed)

#Save
saveRDS(Effect_methanol_depending_length_to_plot, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Effect_methanol_depending_length_to_plot.rds")
#Load
Effect_methanol_depending_length_to_plot<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Effect_methanol_depending_length_to_plot.rds")


ggplot(Effect_methanol_depending_length_to_plot, aes(x=exonic.gene.sizes, y=Levels_of_expression, group=Condition)) +
  geom_line(aes(color=Condition))+
  geom_point(aes(color=Condition))+
  scale_color_manual(values=c('#00BFC4','#F8766D'))+
  theme(legend.position="top")

#Without dots
ggplot(Effect_methanol_depending_length_to_plot, aes(x=exonic.gene.sizes, y=Levels_of_expression, group=Condition)) +
  geom_line(aes(color=Condition))+
  scale_color_manual(values=c('#00BFC4','#F8766D'))+
  theme(legend.position="top")
#Without lines
ggplot(Effect_methanol_depending_length_to_plot, aes(x=exonic.gene.sizes, y=Levels_of_expression, group=Condition)) +
  geom_point(aes(color=Condition))+
  scale_color_manual(values=c('#00BFC4','#F8766D'))+
  geom_smooth(method=lm, aes(color=Condition))+
  theme(legend.position="top")

#-------------------Gene classification in 20 Groups separated by length. Each groups have the equal number of genes-----------------------------------------------------

#I need to order the rows in ascending exonic.gene.sizes
Sorted_Effect_methanol_depending_length_to_plot<-Effect_methanol_depending_length_to_plot[order(Effect_methanol_depending_length_to_plot$exonic.gene.sizes),]

#Create a vector with repeated values grouping the genes in 20 groups (I have each gene 2 times because I have the condition Fixed or Fresh as a column)
Groups<-rep(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),each=1410)
Sorted_Effect_methanol_depending_length_to_plot$Groups<-Groups

#Save
saveRDS(Sorted_Effect_methanol_depending_length_to_plot, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Sorted_Effect_methanol_depending_length_to_plot.rds")
#Load
Sorted_Effect_methanol_depending_length_to_plot<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Sorted_Effect_methanol_depending_length_to_plot.rds")


#Bar plot
ggplot(data=Sorted_Effect_methanol_depending_length_to_plot, aes(x=factor(Groups), y=Levels_of_expression, fill=Condition)) +
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))

#Including a column with log10 from the levels of expression
Sorted_Effect_methanol_depending_length_to_plot$log10<-log10(Sorted_Effect_methanol_depending_length_to_plot$Levels_of_expression)

#Bar plot
ggplot(data=Sorted_Effect_methanol_depending_length_to_plot, aes(x=factor(Groups), y=log10, fill=Condition)) +
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))

#Box plots
P1<-ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=factor(Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot()+labs(title="Levels_of_expression")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))
P2<-ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=factor(Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Levels_of_expression_ without outliers")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+ 
  ylim(0,2)

P3<-ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=factor(Groups), y=log10, fill=Condition)) + 
  geom_boxplot()+labs(title="log10(Levels_of_expression)")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggarrange(P1,P2,P3)



#Box plots with dots
P4<-ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=factor(Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot()+labs(title="Levels_of_expression")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)+ geom_jitter(shape=16, position=position_jitter(0.2))
P5<-ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=factor(Groups), y=log10, fill=Condition)) + 
  geom_boxplot()+labs(title="log10(Levels_of_expression)")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)+ geom_jitter(shape=16, position=position_jitter(0.2))
ggarrange(P4,P5)


#Statistic test comparing Fixed and Fresh in each length group
#To keep the script symple I change the name from the dataframe from "Sorted_Effect_methanol_depending_length_to_plot" to "df"
df<-Sorted_Effect_methanol_depending_length_to_plot
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
B_corrected_p_values<-as.list(p.adjust(p_values,method="bonferroni"))
B_Significative<-list()
for (i in 1:20){B_Significative[[i]]<-
  if (B_corrected_p_values[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values<-as.list(p.adjust(p_values,method="BH"))
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
Length_groups<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")
Length_groups_statistics_summary<-data.frame(Length_groups,p_values,Statistic_W,B_corrected_p_values,B_Significative,BH_corrected_p_values,BH_Significative)
Length_groups_statistics_summary

# Export the summary table
time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Length_groups_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_vs_Fixed", col.names = TRUE, row.names = TRUE, append = FALSE)


#Statistic test comparing the different lengths inside Fresh and fixed
Fresh_and_Fixed<-split(df, df$Condition)
#----------------------------------------------------------------Fresh--------------------------------------------------------------------------------------------

Fresh<-Fresh_and_Fixed[[2]]
Fresh_Exp<-Fresh[,c(-1,-3,-4,-5,-7)]
Fresh_Exp<-split(Fresh_Exp, Fresh$Groups)
Fresh_Exp<-lapply(Fresh_Exp, function(x) x[!(names(x) %in% c("Groups"))])
Fresh_exp_list<-list()
for(i in 1:20) {Fresh_exp_list[[i]]<-unlist(Fresh_Exp[[i]])}

#Statistics
#Mean
Expression_Means_Fresh<-list()
for(i in 1:20) {Expression_Means_Fresh[[i]]<-mean(Fresh_exp_list[[i]], na.rm=TRUE)}
#median
Expression_median_Fresh<-list()
for(i in 1:20) {Expression_median_Fresh[[i]]<-median(Fresh_exp_list[[i]], na.rm=TRUE)}
#standard error of the mean (SEM)
std <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))}
Expression_std_Fresh<-list()
for(i in 1:20) {Expression_std_Fresh[[i]]<-std(Fresh_exp_list[[i]], na.rm=TRUE)}

Expression_Means_Fresh<-paste(Expression_Means_Fresh, collapse = NULL)
Expression_median_Fresh<-paste(Expression_median_Fresh, collapse = NULL)
Expression_std_Fresh<-paste(Expression_std_Fresh, collapse = NULL)

#Results
Length_groups<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")
Fresh_statistics_summary<-data.frame(Length_groups,Expression_Means_Fresh,Expression_median_Fresh,Expression_std_Fresh)
Fresh_statistics_summary

# Export the summary table
write.xlsx(Fresh_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_statistics", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------1----------------------------------------------------------------------------
W_test_t_Fresh_1<-list()
for (i in 1:20){W_test_t_Fresh_1[[i]]<-wilcox.test(Fresh_exp_list[[1]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_1<-list()
for (i in 1:20){p_values_Fresh_1[[i]]<-W_test_t_Fresh_1[[i]][["p.value"]]}
Statistic_W_Fresh_1<-list()
for (i in 1:20){Statistic_W_Fresh_1[[i]]<-W_test_t_Fresh_1[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_1<-as.list(p.adjust(p_values_Fresh_1,method="bonferroni"))
B_Significative_Fresh_1<-list()
for (i in 1:20){B_Significative_Fresh_1[[i]]<-
  if (B_corrected_p_values_Fresh_1[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_1<-as.list(p.adjust(p_values_Fresh_1,method="BH"))
BH_Significative_Fresh_1<-list()
for (i in 1:20){BH_Significative_Fresh_1[[i]]<-
  if (BH_corrected_p_values_Fresh_1[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_1<-paste(p_values_Fresh_1, collapse = NULL)
Statistic_W_Fresh_1<-paste(Statistic_W_Fresh_1, collapse = NULL)
B_corrected_p_values_Fresh_1<-paste(B_corrected_p_values_Fresh_1, collapse = NULL)
B_Significative_Fresh_1<-paste(B_Significative_Fresh_1, collapse = NULL)
BH_corrected_p_values_Fresh_1<-paste(BH_corrected_p_values_Fresh_1, collapse = NULL)
BH_Significative_Fresh_1<-paste(BH_Significative_Fresh_1, collapse = NULL)

#Results
Length_groups<-c("1 vs 1","1 vs 2","1 vs 3","1 vs 4","1 vs 5","1 vs 6","1 vs 7","1 vs 8","1 vs 9","1 vs 10","1 vs 11","1 vs 12","1 vs 13","1 vs 14","1 vs 15","1 vs 16","1 vs 17","1 vs 18","1 vs 19","1 vs 20")
Fresh_group1_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_1,Statistic_W_Fresh_1,B_corrected_p_values_Fresh_1,
                                                         B_Significative_Fresh_1,BH_corrected_p_values_Fresh_1,BH_Significative_Fresh_1)
write.xlsx(Fresh_group1_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_1_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------2----------------------------------------------------------------------------
W_test_t_Fresh_2<-list()
for (i in 1:20){W_test_t_Fresh_2[[i]]<-wilcox.test(Fresh_exp_list[[2]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_2<-list()
for (i in 1:20){p_values_Fresh_2[[i]]<-W_test_t_Fresh_2[[i]][["p.value"]]}
Statistic_W_Fresh_2<-list()
for (i in 1:20){Statistic_W_Fresh_2[[i]]<-W_test_t_Fresh_2[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_2<-as.list(p.adjust(p_values_Fresh_2,method="bonferroni"))
B_Significative_Fresh_2<-list()
for (i in 1:20){B_Significative_Fresh_2[[i]]<-
  if (B_corrected_p_values_Fresh_2[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_2<-as.list(p.adjust(p_values_Fresh_2,method="BH"))
BH_Significative_Fresh_2<-list()
for (i in 1:20){BH_Significative_Fresh_2[[i]]<-
  if (BH_corrected_p_values_Fresh_2[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_2<-paste(p_values_Fresh_2, collapse = NULL)
Statistic_W_Fresh_2<-paste(Statistic_W_Fresh_2, collapse = NULL)
B_corrected_p_values_Fresh_2<-paste(B_corrected_p_values_Fresh_2, collapse = NULL)
B_Significative_Fresh_2<-paste(B_Significative_Fresh_2, collapse = NULL)
BH_corrected_p_values_Fresh_2<-paste(BH_corrected_p_values_Fresh_2, collapse = NULL)
BH_Significative_Fresh_2<-paste(BH_Significative_Fresh_2, collapse = NULL)

#Results
Length_groups<-c("2 vs 1","2 vs 2","2 vs 3","2 vs 4","2 vs 5","2 vs 6","2 vs 7","2 vs 8","2 vs 9","2 vs 10","2 vs 11","2 vs 12","2 vs 13","2 vs 14","2 vs 15","2 vs 16","2 vs 17","2 vs 18","2 vs 19","2 vs 20")
Fresh_group2_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_2,Statistic_W_Fresh_2,B_corrected_p_values_Fresh_2,
                                                         B_Significative_Fresh_2,BH_corrected_p_values_Fresh_2,BH_Significative_Fresh_2)
write.xlsx(Fresh_group2_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_2_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------3----------------------------------------------------------------------------
W_test_t_Fresh_3<-list()
for (i in 1:20){W_test_t_Fresh_3[[i]]<-wilcox.test(Fresh_exp_list[[3]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_3<-list()
for (i in 1:20){p_values_Fresh_3[[i]]<-W_test_t_Fresh_3[[i]][["p.value"]]}
Statistic_W_Fresh_3<-list()
for (i in 1:20){Statistic_W_Fresh_3[[i]]<-W_test_t_Fresh_3[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_3<-as.list(p.adjust(p_values_Fresh_3,method="bonferroni"))
B_Significative_Fresh_3<-list()
for (i in 1:20){B_Significative_Fresh_3[[i]]<-
  if (B_corrected_p_values_Fresh_3[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_3<-as.list(p.adjust(p_values_Fresh_3,method="BH"))
BH_Significative_Fresh_3<-list()
for (i in 1:20){BH_Significative_Fresh_3[[i]]<-
  if (BH_corrected_p_values_Fresh_3[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_3<-paste(p_values_Fresh_3, collapse = NULL)
Statistic_W_Fresh_3<-paste(Statistic_W_Fresh_3, collapse = NULL)
B_corrected_p_values_Fresh_3<-paste(B_corrected_p_values_Fresh_3, collapse = NULL)
B_Significative_Fresh_3<-paste(B_Significative_Fresh_3, collapse = NULL)
BH_corrected_p_values_Fresh_3<-paste(BH_corrected_p_values_Fresh_3, collapse = NULL)
BH_Significative_Fresh_3<-paste(BH_Significative_Fresh_3, collapse = NULL)

#Results
Length_groups<-c("3 vs 1","3 vs 2","3 vs 3","3 vs 4","3 vs 5","3 vs 6","3 vs 7","3 vs 8","3 vs 9","3 vs 10","3 vs 11","3 vs 12","3 vs 13","3 vs 14","3 vs 15","3 vs 16","3 vs 17","3 vs 18","3 vs 19","3 vs 20")
Fresh_group3_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_3,Statistic_W_Fresh_3,B_corrected_p_values_Fresh_3,
                                                         B_Significative_Fresh_3,BH_corrected_p_values_Fresh_3,BH_Significative_Fresh_3)
write.xlsx(Fresh_group3_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_3_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------4----------------------------------------------------------------------------
W_test_t_Fresh_4<-list()
for (i in 1:20){W_test_t_Fresh_4[[i]]<-wilcox.test(Fresh_exp_list[[4]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_4<-list()
for (i in 1:20){p_values_Fresh_4[[i]]<-W_test_t_Fresh_4[[i]][["p.value"]]}
Statistic_W_Fresh_4<-list()
for (i in 1:20){Statistic_W_Fresh_4[[i]]<-W_test_t_Fresh_4[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_4<-as.list(p.adjust(p_values_Fresh_4,method="bonferroni"))
B_Significative_Fresh_4<-list()
for (i in 1:20){B_Significative_Fresh_4[[i]]<-
  if (B_corrected_p_values_Fresh_4[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_4<-as.list(p.adjust(p_values_Fresh_4,method="BH"))
BH_Significative_Fresh_4<-list()
for (i in 1:20){BH_Significative_Fresh_4[[i]]<-
  if (BH_corrected_p_values_Fresh_4[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_4<-paste(p_values_Fresh_4, collapse = NULL)
Statistic_W_Fresh_4<-paste(Statistic_W_Fresh_4, collapse = NULL)
B_corrected_p_values_Fresh_4<-paste(B_corrected_p_values_Fresh_4, collapse = NULL)
B_Significative_Fresh_4<-paste(B_Significative_Fresh_4, collapse = NULL)
BH_corrected_p_values_Fresh_4<-paste(BH_corrected_p_values_Fresh_4, collapse = NULL)
BH_Significative_Fresh_4<-paste(BH_Significative_Fresh_4, collapse = NULL)

#Results
Length_groups<-c("4 vs 1","4 vs 2","4 vs 3","4 vs 4","4 vs 5","4 vs 6","4 vs 7","4 vs 8","4 vs 9","4 vs 10","4 vs 11","4 vs 12","4 vs 13","4 vs 14","4 vs 15","4 vs 16","4 vs 17","4 vs 18","4 vs 19","4 vs 20")
Fresh_group4_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_4,Statistic_W_Fresh_4,B_corrected_p_values_Fresh_4,
                                                         B_Significative_Fresh_4,BH_corrected_p_values_Fresh_4,BH_Significative_Fresh_4)
write.xlsx(Fresh_group4_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_4_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------5----------------------------------------------------------------------------
W_test_t_Fresh_5<-list()
for (i in 1:20){W_test_t_Fresh_5[[i]]<-wilcox.test(Fresh_exp_list[[5]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_5<-list()
for (i in 1:20){p_values_Fresh_5[[i]]<-W_test_t_Fresh_5[[i]][["p.value"]]}
Statistic_W_Fresh_5<-list()
for (i in 1:20){Statistic_W_Fresh_5[[i]]<-W_test_t_Fresh_5[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_5<-as.list(p.adjust(p_values_Fresh_5,method="bonferroni"))
B_Significative_Fresh_5<-list()
for (i in 1:20){B_Significative_Fresh_5[[i]]<-
  if (B_corrected_p_values_Fresh_5[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_5<-as.list(p.adjust(p_values_Fresh_5,method="BH"))
BH_Significative_Fresh_5<-list()
for (i in 1:20){BH_Significative_Fresh_5[[i]]<-
  if (BH_corrected_p_values_Fresh_5[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_5<-paste(p_values_Fresh_5, collapse = NULL)
Statistic_W_Fresh_5<-paste(Statistic_W_Fresh_5, collapse = NULL)
B_corrected_p_values_Fresh_5<-paste(B_corrected_p_values_Fresh_5, collapse = NULL)
B_Significative_Fresh_5<-paste(B_Significative_Fresh_5, collapse = NULL)
BH_corrected_p_values_Fresh_5<-paste(BH_corrected_p_values_Fresh_5, collapse = NULL)
BH_Significative_Fresh_5<-paste(BH_Significative_Fresh_5, collapse = NULL)

#Results
Length_groups<-c("5 vs 1","5 vs 2","5 vs 3","5 vs 4","5 vs 5","5 vs 6","5 vs 7","5 vs 8","5 vs 9","5 vs 10","5 vs 11","5 vs 12","5 vs 13","5 vs 14","5 vs 15","5 vs 16","5 vs 17","5 vs 18","5 vs 19","5 vs 20")
Fresh_group5_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_5,Statistic_W_Fresh_5,B_corrected_p_values_Fresh_5,
                                                         B_Significative_Fresh_5,BH_corrected_p_values_Fresh_5,BH_Significative_Fresh_5)
write.xlsx(Fresh_group5_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_5_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------6----------------------------------------------------------------------------
W_test_t_Fresh_6<-list()
for (i in 1:20){W_test_t_Fresh_6[[i]]<-wilcox.test(Fresh_exp_list[[6]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_6<-list()
for (i in 1:20){p_values_Fresh_6[[i]]<-W_test_t_Fresh_6[[i]][["p.value"]]}
Statistic_W_Fresh_6<-list()
for (i in 1:20){Statistic_W_Fresh_6[[i]]<-W_test_t_Fresh_6[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_6<-as.list(p.adjust(p_values_Fresh_6,method="bonferroni"))
B_Significative_Fresh_6<-list()
for (i in 1:20){B_Significative_Fresh_6[[i]]<-
  if (B_corrected_p_values_Fresh_6[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_6<-as.list(p.adjust(p_values_Fresh_6,method="BH"))
BH_Significative_Fresh_6<-list()
for (i in 1:20){BH_Significative_Fresh_6[[i]]<-
  if (BH_corrected_p_values_Fresh_6[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_6<-paste(p_values_Fresh_6, collapse = NULL)
Statistic_W_Fresh_6<-paste(Statistic_W_Fresh_6, collapse = NULL)
B_corrected_p_values_Fresh_6<-paste(B_corrected_p_values_Fresh_6, collapse = NULL)
B_Significative_Fresh_6<-paste(B_Significative_Fresh_6, collapse = NULL)
BH_corrected_p_values_Fresh_6<-paste(BH_corrected_p_values_Fresh_6, collapse = NULL)
BH_Significative_Fresh_6<-paste(BH_Significative_Fresh_6, collapse = NULL)

#Results
Length_groups<-c("6 vs 1","6 vs 2","6 vs 3","6 vs 4","6 vs 5","6 vs 6","6 vs 7","6 vs 8","6 vs 9","6 vs 10","6 vs 11","6 vs 12","6 vs 13","6 vs 14","6 vs 15","6 vs 16","6 vs 17","6 vs 18","6 vs 19","6 vs 20")
Fresh_group6_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_6,Statistic_W_Fresh_6,B_corrected_p_values_Fresh_6,
                                                         B_Significative_Fresh_6,BH_corrected_p_values_Fresh_6,BH_Significative_Fresh_6)
write.xlsx(Fresh_group6_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_6_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------7----------------------------------------------------------------------------
W_test_t_Fresh_7<-list()
for (i in 1:20){W_test_t_Fresh_7[[i]]<-wilcox.test(Fresh_exp_list[[7]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_7<-list()
for (i in 1:20){p_values_Fresh_7[[i]]<-W_test_t_Fresh_7[[i]][["p.value"]]}
Statistic_W_Fresh_7<-list()
for (i in 1:20){Statistic_W_Fresh_7[[i]]<-W_test_t_Fresh_7[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_7<-as.list(p.adjust(p_values_Fresh_7,method="bonferroni"))
B_Significative_Fresh_7<-list()
for (i in 1:20){B_Significative_Fresh_7[[i]]<-
  if (B_corrected_p_values_Fresh_7[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_7<-as.list(p.adjust(p_values_Fresh_7,method="BH"))
BH_Significative_Fresh_7<-list()
for (i in 1:20){BH_Significative_Fresh_7[[i]]<-
  if (BH_corrected_p_values_Fresh_7[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_7<-paste(p_values_Fresh_7, collapse = NULL)
Statistic_W_Fresh_7<-paste(Statistic_W_Fresh_7, collapse = NULL)
B_corrected_p_values_Fresh_7<-paste(B_corrected_p_values_Fresh_7, collapse = NULL)
B_Significative_Fresh_7<-paste(B_Significative_Fresh_7, collapse = NULL)
BH_corrected_p_values_Fresh_7<-paste(BH_corrected_p_values_Fresh_7, collapse = NULL)
BH_Significative_Fresh_7<-paste(BH_Significative_Fresh_7, collapse = NULL)

#Results
Length_groups<-c("7 vs 1","7 vs 2","7 vs 3","7 vs 4","7 vs 5","7 vs 6","7 vs 7","7 vs 8","7 vs 9","7 vs 10","7 vs 11","7 vs 12","7 vs 13","7 vs 14","7 vs 15","7 vs 16","7 vs 17","7 vs 18","7 vs 19","7 vs 20")
Fresh_group7_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_7,Statistic_W_Fresh_7,B_corrected_p_values_Fresh_7,
                                                         B_Significative_Fresh_7,BH_corrected_p_values_Fresh_7,BH_Significative_Fresh_7)
write.xlsx(Fresh_group7_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_7_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------8----------------------------------------------------------------------------
W_test_t_Fresh_8<-list()
for (i in 1:20){W_test_t_Fresh_8[[i]]<-wilcox.test(Fresh_exp_list[[8]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_8<-list()
for (i in 1:20){p_values_Fresh_8[[i]]<-W_test_t_Fresh_8[[i]][["p.value"]]}
Statistic_W_Fresh_8<-list()
for (i in 1:20){Statistic_W_Fresh_8[[i]]<-W_test_t_Fresh_8[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_8<-as.list(p.adjust(p_values_Fresh_8,method="bonferroni"))
B_Significative_Fresh_8<-list()
for (i in 1:20){B_Significative_Fresh_8[[i]]<-
  if (B_corrected_p_values_Fresh_8[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_8<-as.list(p.adjust(p_values_Fresh_8,method="BH"))
BH_Significative_Fresh_8<-list()
for (i in 1:20){BH_Significative_Fresh_8[[i]]<-
  if (BH_corrected_p_values_Fresh_8[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_8<-paste(p_values_Fresh_8, collapse = NULL)
Statistic_W_Fresh_8<-paste(Statistic_W_Fresh_8, collapse = NULL)
B_corrected_p_values_Fresh_8<-paste(B_corrected_p_values_Fresh_8, collapse = NULL)
B_Significative_Fresh_8<-paste(B_Significative_Fresh_8, collapse = NULL)
BH_corrected_p_values_Fresh_8<-paste(BH_corrected_p_values_Fresh_8, collapse = NULL)
BH_Significative_Fresh_8<-paste(BH_Significative_Fresh_8, collapse = NULL)

#Results
Length_groups<-c("8 vs 1","8 vs 2","8 vs 3","8 vs 4","8 vs 5","8 vs 6","8 vs 7","8 vs 8","8 vs 9","8 vs 10","8 vs 11","8 vs 12","8 vs 13","8 vs 14","8 vs 15","8 vs 16","8 vs 17","8 vs 18","8 vs 19","8 vs 20")
Fresh_group8_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_8,Statistic_W_Fresh_8,B_corrected_p_values_Fresh_8,
                                                         B_Significative_Fresh_8,BH_corrected_p_values_Fresh_8,BH_Significative_Fresh_8)
write.xlsx(Fresh_group8_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_8_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------9----------------------------------------------------------------------------
W_test_t_Fresh_9<-list()
for (i in 1:20){W_test_t_Fresh_9[[i]]<-wilcox.test(Fresh_exp_list[[9]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_9<-list()
for (i in 1:20){p_values_Fresh_9[[i]]<-W_test_t_Fresh_9[[i]][["p.value"]]}
Statistic_W_Fresh_9<-list()
for (i in 1:20){Statistic_W_Fresh_9[[i]]<-W_test_t_Fresh_9[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_9<-as.list(p.adjust(p_values_Fresh_9,method="bonferroni"))
B_Significative_Fresh_9<-list()
for (i in 1:20){B_Significative_Fresh_9[[i]]<-
  if (B_corrected_p_values_Fresh_9[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_9<-as.list(p.adjust(p_values_Fresh_9,method="BH"))
BH_Significative_Fresh_9<-list()
for (i in 1:20){BH_Significative_Fresh_9[[i]]<-
  if (BH_corrected_p_values_Fresh_9[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_9<-paste(p_values_Fresh_9, collapse = NULL)
Statistic_W_Fresh_9<-paste(Statistic_W_Fresh_9, collapse = NULL)
B_corrected_p_values_Fresh_9<-paste(B_corrected_p_values_Fresh_9, collapse = NULL)
B_Significative_Fresh_9<-paste(B_Significative_Fresh_9, collapse = NULL)
BH_corrected_p_values_Fresh_9<-paste(BH_corrected_p_values_Fresh_9, collapse = NULL)
BH_Significative_Fresh_9<-paste(BH_Significative_Fresh_9, collapse = NULL)

#Results
Length_groups<-c("9 vs 1","9 vs 2","9 vs 3","9 vs 4","9 vs 5","9 vs 6","9 vs 7","9 vs 8","9 vs 9","9 vs 10","9 vs 11","9 vs 12","9 vs 13","9 vs 14","9 vs 15","9 vs 16","9 vs 17","9 vs 18","9 vs 19","9 vs 20")
Fresh_group9_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_9,Statistic_W_Fresh_9,B_corrected_p_values_Fresh_9,
                                                         B_Significative_Fresh_9,BH_corrected_p_values_Fresh_9,BH_Significative_Fresh_9)
write.xlsx(Fresh_group9_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_9_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------10----------------------------------------------------------------------------
W_test_t_Fresh_10<-list()
for (i in 1:20){W_test_t_Fresh_10[[i]]<-wilcox.test(Fresh_exp_list[[10]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_10<-list()
for (i in 1:20){p_values_Fresh_10[[i]]<-W_test_t_Fresh_10[[i]][["p.value"]]}
Statistic_W_Fresh_10<-list()
for (i in 1:20){Statistic_W_Fresh_10[[i]]<-W_test_t_Fresh_10[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_10<-as.list(p.adjust(p_values_Fresh_10,method="bonferroni"))
B_Significative_Fresh_10<-list()
for (i in 1:20){B_Significative_Fresh_10[[i]]<-
  if (B_corrected_p_values_Fresh_10[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_10<-as.list(p.adjust(p_values_Fresh_10,method="BH"))
BH_Significative_Fresh_10<-list()
for (i in 1:20){BH_Significative_Fresh_10[[i]]<-
  if (BH_corrected_p_values_Fresh_10[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_10<-paste(p_values_Fresh_10, collapse = NULL)
Statistic_W_Fresh_10<-paste(Statistic_W_Fresh_10, collapse = NULL)
B_corrected_p_values_Fresh_10<-paste(B_corrected_p_values_Fresh_10, collapse = NULL)
B_Significative_Fresh_10<-paste(B_Significative_Fresh_10, collapse = NULL)
BH_corrected_p_values_Fresh_10<-paste(BH_corrected_p_values_Fresh_10, collapse = NULL)
BH_Significative_Fresh_10<-paste(BH_Significative_Fresh_10, collapse = NULL)

#Results
Length_groups<-c("10 vs 1","10 vs 2","10 vs 3","10 vs 4","10 vs 5","10 vs 6","10 vs 7","10 vs 8","10 vs 9","10 vs 10","10 vs 11","10 vs 12","10 vs 13","10 vs 14","10 vs 15","10 vs 16","10 vs 17","10 vs 18","10 vs 19","10 vs 20")
Fresh_group10_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_10,Statistic_W_Fresh_10,B_corrected_p_values_Fresh_10,
                                                         B_Significative_Fresh_10,BH_corrected_p_values_Fresh_10,BH_Significative_Fresh_10)
write.xlsx(Fresh_group10_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_10_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------11----------------------------------------------------------------------------
W_test_t_Fresh_11<-list()
for (i in 1:20){W_test_t_Fresh_11[[i]]<-wilcox.test(Fresh_exp_list[[11]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_11<-list()
for (i in 1:20){p_values_Fresh_11[[i]]<-W_test_t_Fresh_11[[i]][["p.value"]]}
Statistic_W_Fresh_11<-list()
for (i in 1:20){Statistic_W_Fresh_11[[i]]<-W_test_t_Fresh_11[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_11<-as.list(p.adjust(p_values_Fresh_11,method="bonferroni"))
B_Significative_Fresh_11<-list()
for (i in 1:20){B_Significative_Fresh_11[[i]]<-
  if (B_corrected_p_values_Fresh_11[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_11<-as.list(p.adjust(p_values_Fresh_11,method="BH"))
BH_Significative_Fresh_11<-list()
for (i in 1:20){BH_Significative_Fresh_11[[i]]<-
  if (BH_corrected_p_values_Fresh_11[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_11<-paste(p_values_Fresh_11, collapse = NULL)
Statistic_W_Fresh_11<-paste(Statistic_W_Fresh_11, collapse = NULL)
B_corrected_p_values_Fresh_11<-paste(B_corrected_p_values_Fresh_11, collapse = NULL)
B_Significative_Fresh_11<-paste(B_Significative_Fresh_11, collapse = NULL)
BH_corrected_p_values_Fresh_11<-paste(BH_corrected_p_values_Fresh_11, collapse = NULL)
BH_Significative_Fresh_11<-paste(BH_Significative_Fresh_11, collapse = NULL)

#Results
Length_groups<-c("11 vs 1","11 vs 2","11 vs 3","11 vs 4","11 vs 5","11 vs 6","11 vs 7","11 vs 8","11 vs 9","11 vs 10","11 vs 11","11 vs 12","11 vs 13","11 vs 14","11 vs 15","11 vs 16","11 vs 17","11 vs 18","11 vs 19","11 vs 20")
Fresh_group11_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_11,Statistic_W_Fresh_11,B_corrected_p_values_Fresh_11,
                                                         B_Significative_Fresh_11,BH_corrected_p_values_Fresh_11,BH_Significative_Fresh_11)
write.xlsx(Fresh_group11_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_11_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------12----------------------------------------------------------------------------
W_test_t_Fresh_12<-list()
for (i in 1:20){W_test_t_Fresh_12[[i]]<-wilcox.test(Fresh_exp_list[[12]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_12<-list()
for (i in 1:20){p_values_Fresh_12[[i]]<-W_test_t_Fresh_12[[i]][["p.value"]]}
Statistic_W_Fresh_12<-list()
for (i in 1:20){Statistic_W_Fresh_12[[i]]<-W_test_t_Fresh_12[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_12<-as.list(p.adjust(p_values_Fresh_12,method="bonferroni"))
B_Significative_Fresh_12<-list()
for (i in 1:20){B_Significative_Fresh_12[[i]]<-
  if (B_corrected_p_values_Fresh_12[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_12<-as.list(p.adjust(p_values_Fresh_12,method="BH"))
BH_Significative_Fresh_12<-list()
for (i in 1:20){BH_Significative_Fresh_12[[i]]<-
  if (BH_corrected_p_values_Fresh_12[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_12<-paste(p_values_Fresh_12, collapse = NULL)
Statistic_W_Fresh_12<-paste(Statistic_W_Fresh_12, collapse = NULL)
B_corrected_p_values_Fresh_12<-paste(B_corrected_p_values_Fresh_12, collapse = NULL)
B_Significative_Fresh_12<-paste(B_Significative_Fresh_12, collapse = NULL)
BH_corrected_p_values_Fresh_12<-paste(BH_corrected_p_values_Fresh_12, collapse = NULL)
BH_Significative_Fresh_12<-paste(BH_Significative_Fresh_12, collapse = NULL)

#Results
Length_groups<-c("12 vs 1","12 vs 2","12 vs 3","12 vs 4","12 vs 5","12 vs 6","12 vs 7","12 vs 8","12 vs 9","12 vs 10","12 vs 11","12 vs 12","12 vs 13","12 vs 14","12 vs 15","12 vs 16","12 vs 17","12 vs 18","12 vs 19","12 vs 20")
Fresh_group12_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_12,Statistic_W_Fresh_12,B_corrected_p_values_Fresh_12,
                                                          B_Significative_Fresh_12,BH_corrected_p_values_Fresh_12,BH_Significative_Fresh_12)
write.xlsx(Fresh_group12_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_12_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)
#-----------------------------------------------------------------13----------------------------------------------------------------------------
W_test_t_Fresh_13<-list()
for (i in 1:20){W_test_t_Fresh_13[[i]]<-wilcox.test(Fresh_exp_list[[13]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_13<-list()
for (i in 1:20){p_values_Fresh_13[[i]]<-W_test_t_Fresh_13[[i]][["p.value"]]}
Statistic_W_Fresh_13<-list()
for (i in 1:20){Statistic_W_Fresh_13[[i]]<-W_test_t_Fresh_13[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_13<-as.list(p.adjust(p_values_Fresh_13,method="bonferroni"))
B_Significative_Fresh_13<-list()
for (i in 1:20){B_Significative_Fresh_13[[i]]<-
  if (B_corrected_p_values_Fresh_13[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_13<-as.list(p.adjust(p_values_Fresh_13,method="BH"))
BH_Significative_Fresh_13<-list()
for (i in 1:20){BH_Significative_Fresh_13[[i]]<-
  if (BH_corrected_p_values_Fresh_13[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_13<-paste(p_values_Fresh_13, collapse = NULL)
Statistic_W_Fresh_13<-paste(Statistic_W_Fresh_13, collapse = NULL)
B_corrected_p_values_Fresh_13<-paste(B_corrected_p_values_Fresh_13, collapse = NULL)
B_Significative_Fresh_13<-paste(B_Significative_Fresh_13, collapse = NULL)
BH_corrected_p_values_Fresh_13<-paste(BH_corrected_p_values_Fresh_13, collapse = NULL)
BH_Significative_Fresh_13<-paste(BH_Significative_Fresh_13, collapse = NULL)

#Results
Length_groups<-c("13 vs 1","13 vs 2","13 vs 3","13 vs 4","13 vs 5","13 vs 6","13 vs 7","13 vs 8","13 vs 9","13 vs 10","13 vs 11","13 vs 12","13 vs 13","13 vs 14","13 vs 15","13 vs 16","13 vs 17","13 vs 18","13 vs 19","13 vs 20")
Fresh_group13_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_13,Statistic_W_Fresh_13,B_corrected_p_values_Fresh_13,
                                                          B_Significative_Fresh_13,BH_corrected_p_values_Fresh_13,BH_Significative_Fresh_13)
write.xlsx(Fresh_group13_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_13_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)
#-----------------------------------------------------------------14----------------------------------------------------------------------------
W_test_t_Fresh_14<-list()
for (i in 1:20){W_test_t_Fresh_14[[i]]<-wilcox.test(Fresh_exp_list[[14]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_14<-list()
for (i in 1:20){p_values_Fresh_14[[i]]<-W_test_t_Fresh_14[[i]][["p.value"]]}
Statistic_W_Fresh_14<-list()
for (i in 1:20){Statistic_W_Fresh_14[[i]]<-W_test_t_Fresh_14[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_14<-as.list(p.adjust(p_values_Fresh_14,method="bonferroni"))
B_Significative_Fresh_14<-list()
for (i in 1:20){B_Significative_Fresh_14[[i]]<-
  if (B_corrected_p_values_Fresh_14[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_14<-as.list(p.adjust(p_values_Fresh_14,method="BH"))
BH_Significative_Fresh_14<-list()
for (i in 1:20){BH_Significative_Fresh_14[[i]]<-
  if (BH_corrected_p_values_Fresh_14[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_14<-paste(p_values_Fresh_14, collapse = NULL)
Statistic_W_Fresh_14<-paste(Statistic_W_Fresh_14, collapse = NULL)
B_corrected_p_values_Fresh_14<-paste(B_corrected_p_values_Fresh_14, collapse = NULL)
B_Significative_Fresh_14<-paste(B_Significative_Fresh_14, collapse = NULL)
BH_corrected_p_values_Fresh_14<-paste(BH_corrected_p_values_Fresh_14, collapse = NULL)
BH_Significative_Fresh_14<-paste(BH_Significative_Fresh_14, collapse = NULL)

#Results
Length_groups<-c("14 vs 1","14 vs 2","14 vs 3","14 vs 4","14 vs 5","14 vs 6","14 vs 7","14 vs 8","14 vs 9","14 vs 10","14 vs 11","14 vs 12","14 vs 13","14 vs 14","14 vs 15","14 vs 16","14 vs 17","14 vs 18","14 vs 19","14 vs 20")
Fresh_group14_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_14,Statistic_W_Fresh_14,B_corrected_p_values_Fresh_14,
                                                          B_Significative_Fresh_14,BH_corrected_p_values_Fresh_14,BH_Significative_Fresh_14)
write.xlsx(Fresh_group14_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_14_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------15----------------------------------------------------------------------------
W_test_t_Fresh_15<-list()
for (i in 1:20){W_test_t_Fresh_15[[i]]<-wilcox.test(Fresh_exp_list[[15]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_15<-list()
for (i in 1:20){p_values_Fresh_15[[i]]<-W_test_t_Fresh_15[[i]][["p.value"]]}
Statistic_W_Fresh_15<-list()
for (i in 1:20){Statistic_W_Fresh_15[[i]]<-W_test_t_Fresh_15[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_15<-as.list(p.adjust(p_values_Fresh_15,method="bonferroni"))
B_Significative_Fresh_15<-list()
for (i in 1:20){B_Significative_Fresh_15[[i]]<-
  if (B_corrected_p_values_Fresh_15[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_15<-as.list(p.adjust(p_values_Fresh_15,method="BH"))
BH_Significative_Fresh_15<-list()
for (i in 1:20){BH_Significative_Fresh_15[[i]]<-
  if (BH_corrected_p_values_Fresh_15[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_15<-paste(p_values_Fresh_15, collapse = NULL)
Statistic_W_Fresh_15<-paste(Statistic_W_Fresh_15, collapse = NULL)
B_corrected_p_values_Fresh_15<-paste(B_corrected_p_values_Fresh_15, collapse = NULL)
B_Significative_Fresh_15<-paste(B_Significative_Fresh_15, collapse = NULL)
BH_corrected_p_values_Fresh_15<-paste(BH_corrected_p_values_Fresh_15, collapse = NULL)
BH_Significative_Fresh_15<-paste(BH_Significative_Fresh_15, collapse = NULL)

#Results
Length_groups<-c("15 vs 1","15 vs 2","15 vs 3","15 vs 4","15 vs 5","15 vs 6","15 vs 7","15 vs 8","15 vs 9","15 vs 10","15 vs 11","15 vs 12","15 vs 13","15 vs 14","15 vs 15","15 vs 16","15 vs 17","15 vs 18","15 vs 19","15 vs 20")
Fresh_group15_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_15,Statistic_W_Fresh_15,B_corrected_p_values_Fresh_15,
                                                          B_Significative_Fresh_15,BH_corrected_p_values_Fresh_15,BH_Significative_Fresh_15)
write.xlsx(Fresh_group15_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_15_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)
#-----------------------------------------------------------------16----------------------------------------------------------------------------
W_test_t_Fresh_16<-list()
for (i in 1:20){W_test_t_Fresh_16[[i]]<-wilcox.test(Fresh_exp_list[[16]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_16<-list()
for (i in 1:20){p_values_Fresh_16[[i]]<-W_test_t_Fresh_16[[i]][["p.value"]]}
Statistic_W_Fresh_16<-list()
for (i in 1:20){Statistic_W_Fresh_16[[i]]<-W_test_t_Fresh_16[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_16<-as.list(p.adjust(p_values_Fresh_16,method="bonferroni"))
B_Significative_Fresh_16<-list()
for (i in 1:20){B_Significative_Fresh_16[[i]]<-
  if (B_corrected_p_values_Fresh_16[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_16<-as.list(p.adjust(p_values_Fresh_16,method="BH"))
BH_Significative_Fresh_16<-list()
for (i in 1:20){BH_Significative_Fresh_16[[i]]<-
  if (BH_corrected_p_values_Fresh_16[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_16<-paste(p_values_Fresh_16, collapse = NULL)
Statistic_W_Fresh_16<-paste(Statistic_W_Fresh_16, collapse = NULL)
B_corrected_p_values_Fresh_16<-paste(B_corrected_p_values_Fresh_16, collapse = NULL)
B_Significative_Fresh_16<-paste(B_Significative_Fresh_16, collapse = NULL)
BH_corrected_p_values_Fresh_16<-paste(BH_corrected_p_values_Fresh_16, collapse = NULL)
BH_Significative_Fresh_16<-paste(BH_Significative_Fresh_16, collapse = NULL)

#Results
Length_groups<-c("16 vs 1","16 vs 2","16 vs 3","16 vs 4","16 vs 5","16 vs 6","16 vs 7","16 vs 8","16 vs 9","16 vs 10","16 vs 11","16 vs 12","16 vs 13","16 vs 14","16 vs 15","16 vs 16","16 vs 17","16 vs 18","16 vs 19","16 vs 20")
Fresh_group16_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_16,Statistic_W_Fresh_16,B_corrected_p_values_Fresh_16,
                                                          B_Significative_Fresh_16,BH_corrected_p_values_Fresh_16,BH_Significative_Fresh_16)
write.xlsx(Fresh_group16_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_16_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------17----------------------------------------------------------------------------
W_test_t_Fresh_17<-list()
for (i in 1:20){W_test_t_Fresh_17[[i]]<-wilcox.test(Fresh_exp_list[[17]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_17<-list()
for (i in 1:20){p_values_Fresh_17[[i]]<-W_test_t_Fresh_17[[i]][["p.value"]]}
Statistic_W_Fresh_17<-list()
for (i in 1:20){Statistic_W_Fresh_17[[i]]<-W_test_t_Fresh_17[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_17<-as.list(p.adjust(p_values_Fresh_17,method="bonferroni"))
B_Significative_Fresh_17<-list()
for (i in 1:20){B_Significative_Fresh_17[[i]]<-
  if (B_corrected_p_values_Fresh_17[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_17<-as.list(p.adjust(p_values_Fresh_17,method="BH"))
BH_Significative_Fresh_17<-list()
for (i in 1:20){BH_Significative_Fresh_17[[i]]<-
  if (BH_corrected_p_values_Fresh_17[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_17<-paste(p_values_Fresh_17, collapse = NULL)
Statistic_W_Fresh_17<-paste(Statistic_W_Fresh_17, collapse = NULL)
B_corrected_p_values_Fresh_17<-paste(B_corrected_p_values_Fresh_17, collapse = NULL)
B_Significative_Fresh_17<-paste(B_Significative_Fresh_17, collapse = NULL)
BH_corrected_p_values_Fresh_17<-paste(BH_corrected_p_values_Fresh_17, collapse = NULL)
BH_Significative_Fresh_17<-paste(BH_Significative_Fresh_17, collapse = NULL)

#Results
Length_groups<-c("17 vs 1","17 vs 2","17 vs 3","17 vs 4","17 vs 5","17 vs 6","17 vs 7","17 vs 8","17 vs 9","17 vs 10","17 vs 11","17 vs 12","17 vs 13","17 vs 14","17 vs 15","17 vs 16","17 vs 17","17 vs 18","17 vs 19","17 vs 20")
Fresh_group17_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_17,Statistic_W_Fresh_17,B_corrected_p_values_Fresh_17,
                                                          B_Significative_Fresh_17,BH_corrected_p_values_Fresh_17,BH_Significative_Fresh_17)
write.xlsx(Fresh_group17_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_17_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)


#-----------------------------------------------------------------18----------------------------------------------------------------------------
W_test_t_Fresh_18<-list()
for (i in 1:20){W_test_t_Fresh_18[[i]]<-wilcox.test(Fresh_exp_list[[18]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_18<-list()
for (i in 1:20){p_values_Fresh_18[[i]]<-W_test_t_Fresh_18[[i]][["p.value"]]}
Statistic_W_Fresh_18<-list()
for (i in 1:20){Statistic_W_Fresh_18[[i]]<-W_test_t_Fresh_18[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_18<-as.list(p.adjust(p_values_Fresh_18,method="bonferroni"))
B_Significative_Fresh_18<-list()
for (i in 1:20){B_Significative_Fresh_18[[i]]<-
  if (B_corrected_p_values_Fresh_18[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_18<-as.list(p.adjust(p_values_Fresh_18,method="BH"))
BH_Significative_Fresh_18<-list()
for (i in 1:20){BH_Significative_Fresh_18[[i]]<-
  if (BH_corrected_p_values_Fresh_18[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_18<-paste(p_values_Fresh_18, collapse = NULL)
Statistic_W_Fresh_18<-paste(Statistic_W_Fresh_18, collapse = NULL)
B_corrected_p_values_Fresh_18<-paste(B_corrected_p_values_Fresh_18, collapse = NULL)
B_Significative_Fresh_18<-paste(B_Significative_Fresh_18, collapse = NULL)
BH_corrected_p_values_Fresh_18<-paste(BH_corrected_p_values_Fresh_18, collapse = NULL)
BH_Significative_Fresh_18<-paste(BH_Significative_Fresh_18, collapse = NULL)

#Results
Length_groups<-c("18 vs 1","18 vs 2","18 vs 3","18 vs 4","18 vs 5","18 vs 6","18 vs 7","18 vs 8","18 vs 9","18 vs 10","18 vs 11","18 vs 12","18 vs 13","18 vs 14","18 vs 15","18 vs 16","18 vs 17","18 vs 18","18 vs 19","18 vs 20")
Fresh_group18_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_18,Statistic_W_Fresh_18,B_corrected_p_values_Fresh_18,
                                                          B_Significative_Fresh_18,BH_corrected_p_values_Fresh_18,BH_Significative_Fresh_18)
write.xlsx(Fresh_group18_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_18_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------19----------------------------------------------------------------------------
W_test_t_Fresh_19<-list()
for (i in 1:20){W_test_t_Fresh_19[[i]]<-wilcox.test(Fresh_exp_list[[19]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_19<-list()
for (i in 1:20){p_values_Fresh_19[[i]]<-W_test_t_Fresh_19[[i]][["p.value"]]}
Statistic_W_Fresh_19<-list()
for (i in 1:20){Statistic_W_Fresh_19[[i]]<-W_test_t_Fresh_19[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_19<-as.list(p.adjust(p_values_Fresh_19,method="bonferroni"))
B_Significative_Fresh_19<-list()
for (i in 1:20){B_Significative_Fresh_19[[i]]<-
  if (B_corrected_p_values_Fresh_19[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_19<-as.list(p.adjust(p_values_Fresh_19,method="BH"))
BH_Significative_Fresh_19<-list()
for (i in 1:20){BH_Significative_Fresh_19[[i]]<-
  if (BH_corrected_p_values_Fresh_19[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_19<-paste(p_values_Fresh_19, collapse = NULL)
Statistic_W_Fresh_19<-paste(Statistic_W_Fresh_19, collapse = NULL)
B_corrected_p_values_Fresh_19<-paste(B_corrected_p_values_Fresh_19, collapse = NULL)
B_Significative_Fresh_19<-paste(B_Significative_Fresh_19, collapse = NULL)
BH_corrected_p_values_Fresh_19<-paste(BH_corrected_p_values_Fresh_19, collapse = NULL)
BH_Significative_Fresh_19<-paste(BH_Significative_Fresh_19, collapse = NULL)

#Results
Length_groups<-c("19 vs 1","19 vs 2","19 vs 3","19 vs 4","19 vs 5","19 vs 6","19 vs 7","19 vs 8","19 vs 9","19 vs 10","19 vs 11","19 vs 12","19 vs 13","19 vs 14","19 vs 15","19 vs 16","19 vs 17","19 vs 18","19 vs 19","19 vs 20")
Fresh_group19_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_19,Statistic_W_Fresh_19,B_corrected_p_values_Fresh_19,
                                                          B_Significative_Fresh_19,BH_corrected_p_values_Fresh_19,BH_Significative_Fresh_19)
write.xlsx(Fresh_group19_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_19_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)
#-----------------------------------------------------------------20----------------------------------------------------------------------------
W_test_t_Fresh_20<-list()
for (i in 1:20){W_test_t_Fresh_20[[i]]<-wilcox.test(Fresh_exp_list[[20]],Fresh_exp_list[[i]], paired=FALSE)}
p_values_Fresh_20<-list()
for (i in 1:20){p_values_Fresh_20[[i]]<-W_test_t_Fresh_20[[i]][["p.value"]]}
Statistic_W_Fresh_20<-list()
for (i in 1:20){Statistic_W_Fresh_20[[i]]<-W_test_t_Fresh_20[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fresh_20<-as.list(p.adjust(p_values_Fresh_20,method="bonferroni"))
B_Significative_Fresh_20<-list()
for (i in 1:20){B_Significative_Fresh_20[[i]]<-
  if (B_corrected_p_values_Fresh_20[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fresh_20<-as.list(p.adjust(p_values_Fresh_20,method="BH"))
BH_Significative_Fresh_20<-list()
for (i in 1:20){BH_Significative_Fresh_20[[i]]<-
  if (BH_corrected_p_values_Fresh_20[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fresh_20<-paste(p_values_Fresh_20, collapse = NULL)
Statistic_W_Fresh_20<-paste(Statistic_W_Fresh_20, collapse = NULL)
B_corrected_p_values_Fresh_20<-paste(B_corrected_p_values_Fresh_20, collapse = NULL)
B_Significative_Fresh_20<-paste(B_Significative_Fresh_20, collapse = NULL)
BH_corrected_p_values_Fresh_20<-paste(BH_corrected_p_values_Fresh_20, collapse = NULL)
BH_Significative_Fresh_20<-paste(BH_Significative_Fresh_20, collapse = NULL)

#Results
Length_groups<-c("20 vs 1","20 vs 2","20 vs 3","20 vs 4","20 vs 5","20 vs 6","20 vs 7","20 vs 8","20 vs 9","20 vs 10","20 vs 11","20 vs 12","20 vs 13","20 vs 14","20 vs 15","20 vs 16","20 vs 17","20 vs 18","20 vs 19","20 vs 20")
Fresh_group20_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fresh_20,Statistic_W_Fresh_20,B_corrected_p_values_Fresh_20,
                                                          B_Significative_Fresh_20,BH_corrected_p_values_Fresh_20,BH_Significative_Fresh_20)
write.xlsx(Fresh_group20_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fresh_20_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#----------------------------------------------------------------Fixed--------------------------------------------------------------------------------------------

Fixed<-Fresh_and_Fixed[[1]]
Fixed_Exp<-Fixed[,c(-1,-3,-4,-5,-7)]
Fixed_Exp<-split(Fixed_Exp, Fixed$Groups)
Fixed_Exp<-lapply(Fixed_Exp, function(x) x[!(names(x) %in% c("Groups"))])
Fixed_exp_list<-list()
for(i in 1:20) {Fixed_exp_list[[i]]<-unlist(Fixed_Exp[[i]])}

#Statistics
#Mean
Expression_Means_Fixed<-list()
for(i in 1:20) {Expression_Means_Fixed[[i]]<-mean(Fixed_exp_list[[i]], na.rm=TRUE)}
#median
Expression_median_Fixed<-list()
for(i in 1:20) {Expression_median_Fixed[[i]]<-median(Fixed_exp_list[[i]], na.rm=TRUE)}
#standard error of the mean (SEM)
std <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))}
Expression_std_Fixed<-list()
for(i in 1:20) {Expression_std_Fixed[[i]]<-std(Fixed_exp_list[[i]], na.rm=TRUE)}

Expression_Means_Fixed<-paste(Expression_Means_Fixed, collapse = NULL)
Expression_median_Fixed<-paste(Expression_median_Fixed, collapse = NULL)
Expression_std_Fixed<-paste(Expression_std_Fixed, collapse = NULL)

#Results
Length_groups<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")
Fixed_statistics_summary<-data.frame(Length_groups,Expression_Means_Fixed,Expression_median_Fixed,Expression_std_Fixed)
Fixed_statistics_summary

# Export the summary table
write.xlsx(Fixed_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_statistics", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------1----------------------------------------------------------------------------
W_test_t_Fixed_1<-list()
for (i in 1:20){W_test_t_Fixed_1[[i]]<-wilcox.test(Fixed_exp_list[[1]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_1<-list()
for (i in 1:20){p_values_Fixed_1[[i]]<-W_test_t_Fixed_1[[i]][["p.value"]]}
Statistic_W_Fixed_1<-list()
for (i in 1:20){Statistic_W_Fixed_1[[i]]<-W_test_t_Fixed_1[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_1<-as.list(p.adjust(p_values_Fixed_1,method="bonferroni"))
B_Significative_Fixed_1<-list()
for (i in 1:20){B_Significative_Fixed_1[[i]]<-
  if (B_corrected_p_values_Fixed_1[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_1<-as.list(p.adjust(p_values_Fixed_1,method="BH"))
BH_Significative_Fixed_1<-list()
for (i in 1:20){BH_Significative_Fixed_1[[i]]<-
  if (BH_corrected_p_values_Fixed_1[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_1<-paste(p_values_Fixed_1, collapse = NULL)
Statistic_W_Fixed_1<-paste(Statistic_W_Fixed_1, collapse = NULL)
B_corrected_p_values_Fixed_1<-paste(B_corrected_p_values_Fixed_1, collapse = NULL)
B_Significative_Fixed_1<-paste(B_Significative_Fixed_1, collapse = NULL)
BH_corrected_p_values_Fixed_1<-paste(BH_corrected_p_values_Fixed_1, collapse = NULL)
BH_Significative_Fixed_1<-paste(BH_Significative_Fixed_1, collapse = NULL)

#Results
Length_groups<-c("1 vs 1","1 vs 2","1 vs 3","1 vs 4","1 vs 5","1 vs 6","1 vs 7","1 vs 8","1 vs 9","1 vs 10","1 vs 11","1 vs 12","1 vs 13","1 vs 14","1 vs 15","1 vs 16","1 vs 17","1 vs 18","1 vs 19","1 vs 20")
Fixed_group1_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_1,Statistic_W_Fixed_1,B_corrected_p_values_Fixed_1,
                                                         B_Significative_Fixed_1,BH_corrected_p_values_Fixed_1,BH_Significative_Fixed_1)
write.xlsx(Fixed_group1_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_1_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------2----------------------------------------------------------------------------
W_test_t_Fixed_2<-list()
for (i in 1:20){W_test_t_Fixed_2[[i]]<-wilcox.test(Fixed_exp_list[[2]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_2<-list()
for (i in 1:20){p_values_Fixed_2[[i]]<-W_test_t_Fixed_2[[i]][["p.value"]]}
Statistic_W_Fixed_2<-list()
for (i in 1:20){Statistic_W_Fixed_2[[i]]<-W_test_t_Fixed_2[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_2<-as.list(p.adjust(p_values_Fixed_2,method="bonferroni"))
B_Significative_Fixed_2<-list()
for (i in 1:20){B_Significative_Fixed_2[[i]]<-
  if (B_corrected_p_values_Fixed_2[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_2<-as.list(p.adjust(p_values_Fixed_2,method="BH"))
BH_Significative_Fixed_2<-list()
for (i in 1:20){BH_Significative_Fixed_2[[i]]<-
  if (BH_corrected_p_values_Fixed_2[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_2<-paste(p_values_Fixed_2, collapse = NULL)
Statistic_W_Fixed_2<-paste(Statistic_W_Fixed_2, collapse = NULL)
B_corrected_p_values_Fixed_2<-paste(B_corrected_p_values_Fixed_2, collapse = NULL)
B_Significative_Fixed_2<-paste(B_Significative_Fixed_2, collapse = NULL)
BH_corrected_p_values_Fixed_2<-paste(BH_corrected_p_values_Fixed_2, collapse = NULL)
BH_Significative_Fixed_2<-paste(BH_Significative_Fixed_2, collapse = NULL)

#Results
Length_groups<-c("2 vs 1","2 vs 2","2 vs 3","2 vs 4","2 vs 5","2 vs 6","2 vs 7","2 vs 8","2 vs 9","2 vs 10","2 vs 11","2 vs 12","2 vs 13","2 vs 14","2 vs 15","2 vs 16","2 vs 17","2 vs 18","2 vs 19","2 vs 20")
Fixed_group2_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_2,Statistic_W_Fixed_2,B_corrected_p_values_Fixed_2,
                                                         B_Significative_Fixed_2,BH_corrected_p_values_Fixed_2,BH_Significative_Fixed_2)
write.xlsx(Fixed_group2_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_2_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------3----------------------------------------------------------------------------
W_test_t_Fixed_3<-list()
for (i in 1:20){W_test_t_Fixed_3[[i]]<-wilcox.test(Fixed_exp_list[[3]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_3<-list()
for (i in 1:20){p_values_Fixed_3[[i]]<-W_test_t_Fixed_3[[i]][["p.value"]]}
Statistic_W_Fixed_3<-list()
for (i in 1:20){Statistic_W_Fixed_3[[i]]<-W_test_t_Fixed_3[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_3<-as.list(p.adjust(p_values_Fixed_3,method="bonferroni"))
B_Significative_Fixed_3<-list()
for (i in 1:20){B_Significative_Fixed_3[[i]]<-
  if (B_corrected_p_values_Fixed_3[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_3<-as.list(p.adjust(p_values_Fixed_3,method="BH"))
BH_Significative_Fixed_3<-list()
for (i in 1:20){BH_Significative_Fixed_3[[i]]<-
  if (BH_corrected_p_values_Fixed_3[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_3<-paste(p_values_Fixed_3, collapse = NULL)
Statistic_W_Fixed_3<-paste(Statistic_W_Fixed_3, collapse = NULL)
B_corrected_p_values_Fixed_3<-paste(B_corrected_p_values_Fixed_3, collapse = NULL)
B_Significative_Fixed_3<-paste(B_Significative_Fixed_3, collapse = NULL)
BH_corrected_p_values_Fixed_3<-paste(BH_corrected_p_values_Fixed_3, collapse = NULL)
BH_Significative_Fixed_3<-paste(BH_Significative_Fixed_3, collapse = NULL)

#Results
Length_groups<-c("3 vs 1","3 vs 2","3 vs 3","3 vs 4","3 vs 5","3 vs 6","3 vs 7","3 vs 8","3 vs 9","3 vs 10","3 vs 11","3 vs 12","3 vs 13","3 vs 14","3 vs 15","3 vs 16","3 vs 17","3 vs 18","3 vs 19","3 vs 20")
Fixed_group3_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_3,Statistic_W_Fixed_3,B_corrected_p_values_Fixed_3,
                                                         B_Significative_Fixed_3,BH_corrected_p_values_Fixed_3,BH_Significative_Fixed_3)
write.xlsx(Fixed_group3_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_3_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------4----------------------------------------------------------------------------
W_test_t_Fixed_4<-list()
for (i in 1:20){W_test_t_Fixed_4[[i]]<-wilcox.test(Fixed_exp_list[[4]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_4<-list()
for (i in 1:20){p_values_Fixed_4[[i]]<-W_test_t_Fixed_4[[i]][["p.value"]]}
Statistic_W_Fixed_4<-list()
for (i in 1:20){Statistic_W_Fixed_4[[i]]<-W_test_t_Fixed_4[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_4<-as.list(p.adjust(p_values_Fixed_4,method="bonferroni"))
B_Significative_Fixed_4<-list()
for (i in 1:20){B_Significative_Fixed_4[[i]]<-
  if (B_corrected_p_values_Fixed_4[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_4<-as.list(p.adjust(p_values_Fixed_4,method="BH"))
BH_Significative_Fixed_4<-list()
for (i in 1:20){BH_Significative_Fixed_4[[i]]<-
  if (BH_corrected_p_values_Fixed_4[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_4<-paste(p_values_Fixed_4, collapse = NULL)
Statistic_W_Fixed_4<-paste(Statistic_W_Fixed_4, collapse = NULL)
B_corrected_p_values_Fixed_4<-paste(B_corrected_p_values_Fixed_4, collapse = NULL)
B_Significative_Fixed_4<-paste(B_Significative_Fixed_4, collapse = NULL)
BH_corrected_p_values_Fixed_4<-paste(BH_corrected_p_values_Fixed_4, collapse = NULL)
BH_Significative_Fixed_4<-paste(BH_Significative_Fixed_4, collapse = NULL)

#Results
Length_groups<-c("4 vs 1","4 vs 2","4 vs 3","4 vs 4","4 vs 5","4 vs 6","4 vs 7","4 vs 8","4 vs 9","4 vs 10","4 vs 11","4 vs 12","4 vs 13","4 vs 14","4 vs 15","4 vs 16","4 vs 17","4 vs 18","4 vs 19","4 vs 20")
Fixed_group4_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_4,Statistic_W_Fixed_4,B_corrected_p_values_Fixed_4,
                                                         B_Significative_Fixed_4,BH_corrected_p_values_Fixed_4,BH_Significative_Fixed_4)
write.xlsx(Fixed_group4_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_4_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------5----------------------------------------------------------------------------
W_test_t_Fixed_5<-list()
for (i in 1:20){W_test_t_Fixed_5[[i]]<-wilcox.test(Fixed_exp_list[[5]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_5<-list()
for (i in 1:20){p_values_Fixed_5[[i]]<-W_test_t_Fixed_5[[i]][["p.value"]]}
Statistic_W_Fixed_5<-list()
for (i in 1:20){Statistic_W_Fixed_5[[i]]<-W_test_t_Fixed_5[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_5<-as.list(p.adjust(p_values_Fixed_5,method="bonferroni"))
B_Significative_Fixed_5<-list()
for (i in 1:20){B_Significative_Fixed_5[[i]]<-
  if (B_corrected_p_values_Fixed_5[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_5<-as.list(p.adjust(p_values_Fixed_5,method="BH"))
BH_Significative_Fixed_5<-list()
for (i in 1:20){BH_Significative_Fixed_5[[i]]<-
  if (BH_corrected_p_values_Fixed_5[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_5<-paste(p_values_Fixed_5, collapse = NULL)
Statistic_W_Fixed_5<-paste(Statistic_W_Fixed_5, collapse = NULL)
B_corrected_p_values_Fixed_5<-paste(B_corrected_p_values_Fixed_5, collapse = NULL)
B_Significative_Fixed_5<-paste(B_Significative_Fixed_5, collapse = NULL)
BH_corrected_p_values_Fixed_5<-paste(BH_corrected_p_values_Fixed_5, collapse = NULL)
BH_Significative_Fixed_5<-paste(BH_Significative_Fixed_5, collapse = NULL)

#Results
Length_groups<-c("5 vs 1","5 vs 2","5 vs 3","5 vs 4","5 vs 5","5 vs 6","5 vs 7","5 vs 8","5 vs 9","5 vs 10","5 vs 11","5 vs 12","5 vs 13","5 vs 14","5 vs 15","5 vs 16","5 vs 17","5 vs 18","5 vs 19","5 vs 20")
Fixed_group5_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_5,Statistic_W_Fixed_5,B_corrected_p_values_Fixed_5,
                                                         B_Significative_Fixed_5,BH_corrected_p_values_Fixed_5,BH_Significative_Fixed_5)
write.xlsx(Fixed_group5_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_5_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------6----------------------------------------------------------------------------
W_test_t_Fixed_6<-list()
for (i in 1:20){W_test_t_Fixed_6[[i]]<-wilcox.test(Fixed_exp_list[[6]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_6<-list()
for (i in 1:20){p_values_Fixed_6[[i]]<-W_test_t_Fixed_6[[i]][["p.value"]]}
Statistic_W_Fixed_6<-list()
for (i in 1:20){Statistic_W_Fixed_6[[i]]<-W_test_t_Fixed_6[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_6<-as.list(p.adjust(p_values_Fixed_6,method="bonferroni"))
B_Significative_Fixed_6<-list()
for (i in 1:20){B_Significative_Fixed_6[[i]]<-
  if (B_corrected_p_values_Fixed_6[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_6<-as.list(p.adjust(p_values_Fixed_6,method="BH"))
BH_Significative_Fixed_6<-list()
for (i in 1:20){BH_Significative_Fixed_6[[i]]<-
  if (BH_corrected_p_values_Fixed_6[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_6<-paste(p_values_Fixed_6, collapse = NULL)
Statistic_W_Fixed_6<-paste(Statistic_W_Fixed_6, collapse = NULL)
B_corrected_p_values_Fixed_6<-paste(B_corrected_p_values_Fixed_6, collapse = NULL)
B_Significative_Fixed_6<-paste(B_Significative_Fixed_6, collapse = NULL)
BH_corrected_p_values_Fixed_6<-paste(BH_corrected_p_values_Fixed_6, collapse = NULL)
BH_Significative_Fixed_6<-paste(BH_Significative_Fixed_6, collapse = NULL)

#Results
Length_groups<-c("6 vs 1","6 vs 2","6 vs 3","6 vs 4","6 vs 5","6 vs 6","6 vs 7","6 vs 8","6 vs 9","6 vs 10","6 vs 11","6 vs 12","6 vs 13","6 vs 14","6 vs 15","6 vs 16","6 vs 17","6 vs 18","6 vs 19","6 vs 20")
Fixed_group6_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_6,Statistic_W_Fixed_6,B_corrected_p_values_Fixed_6,
                                                         B_Significative_Fixed_6,BH_corrected_p_values_Fixed_6,BH_Significative_Fixed_6)
write.xlsx(Fixed_group6_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_6_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------7----------------------------------------------------------------------------
W_test_t_Fixed_7<-list()
for (i in 1:20){W_test_t_Fixed_7[[i]]<-wilcox.test(Fixed_exp_list[[7]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_7<-list()
for (i in 1:20){p_values_Fixed_7[[i]]<-W_test_t_Fixed_7[[i]][["p.value"]]}
Statistic_W_Fixed_7<-list()
for (i in 1:20){Statistic_W_Fixed_7[[i]]<-W_test_t_Fixed_7[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_7<-as.list(p.adjust(p_values_Fixed_7,method="bonferroni"))
B_Significative_Fixed_7<-list()
for (i in 1:20){B_Significative_Fixed_7[[i]]<-
  if (B_corrected_p_values_Fixed_7[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_7<-as.list(p.adjust(p_values_Fixed_7,method="BH"))
BH_Significative_Fixed_7<-list()
for (i in 1:20){BH_Significative_Fixed_7[[i]]<-
  if (BH_corrected_p_values_Fixed_7[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_7<-paste(p_values_Fixed_7, collapse = NULL)
Statistic_W_Fixed_7<-paste(Statistic_W_Fixed_7, collapse = NULL)
B_corrected_p_values_Fixed_7<-paste(B_corrected_p_values_Fixed_7, collapse = NULL)
B_Significative_Fixed_7<-paste(B_Significative_Fixed_7, collapse = NULL)
BH_corrected_p_values_Fixed_7<-paste(BH_corrected_p_values_Fixed_7, collapse = NULL)
BH_Significative_Fixed_7<-paste(BH_Significative_Fixed_7, collapse = NULL)

#Results
Length_groups<-c("7 vs 1","7 vs 2","7 vs 3","7 vs 4","7 vs 5","7 vs 6","7 vs 7","7 vs 8","7 vs 9","7 vs 10","7 vs 11","7 vs 12","7 vs 13","7 vs 14","7 vs 15","7 vs 16","7 vs 17","7 vs 18","7 vs 19","7 vs 20")
Fixed_group7_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_7,Statistic_W_Fixed_7,B_corrected_p_values_Fixed_7,
                                                         B_Significative_Fixed_7,BH_corrected_p_values_Fixed_7,BH_Significative_Fixed_7)
write.xlsx(Fixed_group7_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_7_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------8----------------------------------------------------------------------------
W_test_t_Fixed_8<-list()
for (i in 1:20){W_test_t_Fixed_8[[i]]<-wilcox.test(Fixed_exp_list[[8]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_8<-list()
for (i in 1:20){p_values_Fixed_8[[i]]<-W_test_t_Fixed_8[[i]][["p.value"]]}
Statistic_W_Fixed_8<-list()
for (i in 1:20){Statistic_W_Fixed_8[[i]]<-W_test_t_Fixed_8[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_8<-as.list(p.adjust(p_values_Fixed_8,method="bonferroni"))
B_Significative_Fixed_8<-list()
for (i in 1:20){B_Significative_Fixed_8[[i]]<-
  if (B_corrected_p_values_Fixed_8[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_8<-as.list(p.adjust(p_values_Fixed_8,method="BH"))
BH_Significative_Fixed_8<-list()
for (i in 1:20){BH_Significative_Fixed_8[[i]]<-
  if (BH_corrected_p_values_Fixed_8[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_8<-paste(p_values_Fixed_8, collapse = NULL)
Statistic_W_Fixed_8<-paste(Statistic_W_Fixed_8, collapse = NULL)
B_corrected_p_values_Fixed_8<-paste(B_corrected_p_values_Fixed_8, collapse = NULL)
B_Significative_Fixed_8<-paste(B_Significative_Fixed_8, collapse = NULL)
BH_corrected_p_values_Fixed_8<-paste(BH_corrected_p_values_Fixed_8, collapse = NULL)
BH_Significative_Fixed_8<-paste(BH_Significative_Fixed_8, collapse = NULL)

#Results
Length_groups<-c("8 vs 1","8 vs 2","8 vs 3","8 vs 4","8 vs 5","8 vs 6","8 vs 7","8 vs 8","8 vs 9","8 vs 10","8 vs 11","8 vs 12","8 vs 13","8 vs 14","8 vs 15","8 vs 16","8 vs 17","8 vs 18","8 vs 19","8 vs 20")
Fixed_group8_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_8,Statistic_W_Fixed_8,B_corrected_p_values_Fixed_8,
                                                         B_Significative_Fixed_8,BH_corrected_p_values_Fixed_8,BH_Significative_Fixed_8)
write.xlsx(Fixed_group8_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_8_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------9----------------------------------------------------------------------------
W_test_t_Fixed_9<-list()
for (i in 1:20){W_test_t_Fixed_9[[i]]<-wilcox.test(Fixed_exp_list[[9]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_9<-list()
for (i in 1:20){p_values_Fixed_9[[i]]<-W_test_t_Fixed_9[[i]][["p.value"]]}
Statistic_W_Fixed_9<-list()
for (i in 1:20){Statistic_W_Fixed_9[[i]]<-W_test_t_Fixed_9[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_9<-as.list(p.adjust(p_values_Fixed_9,method="bonferroni"))
B_Significative_Fixed_9<-list()
for (i in 1:20){B_Significative_Fixed_9[[i]]<-
  if (B_corrected_p_values_Fixed_9[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_9<-as.list(p.adjust(p_values_Fixed_9,method="BH"))
BH_Significative_Fixed_9<-list()
for (i in 1:20){BH_Significative_Fixed_9[[i]]<-
  if (BH_corrected_p_values_Fixed_9[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_9<-paste(p_values_Fixed_9, collapse = NULL)
Statistic_W_Fixed_9<-paste(Statistic_W_Fixed_9, collapse = NULL)
B_corrected_p_values_Fixed_9<-paste(B_corrected_p_values_Fixed_9, collapse = NULL)
B_Significative_Fixed_9<-paste(B_Significative_Fixed_9, collapse = NULL)
BH_corrected_p_values_Fixed_9<-paste(BH_corrected_p_values_Fixed_9, collapse = NULL)
BH_Significative_Fixed_9<-paste(BH_Significative_Fixed_9, collapse = NULL)

#Results
Length_groups<-c("9 vs 1","9 vs 2","9 vs 3","9 vs 4","9 vs 5","9 vs 6","9 vs 7","9 vs 8","9 vs 9","9 vs 10","9 vs 11","9 vs 12","9 vs 13","9 vs 14","9 vs 15","9 vs 16","9 vs 17","9 vs 18","9 vs 19","9 vs 20")
Fixed_group9_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_9,Statistic_W_Fixed_9,B_corrected_p_values_Fixed_9,
                                                         B_Significative_Fixed_9,BH_corrected_p_values_Fixed_9,BH_Significative_Fixed_9)
write.xlsx(Fixed_group9_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_9_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------10----------------------------------------------------------------------------
W_test_t_Fixed_10<-list()
for (i in 1:20){W_test_t_Fixed_10[[i]]<-wilcox.test(Fixed_exp_list[[10]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_10<-list()
for (i in 1:20){p_values_Fixed_10[[i]]<-W_test_t_Fixed_10[[i]][["p.value"]]}
Statistic_W_Fixed_10<-list()
for (i in 1:20){Statistic_W_Fixed_10[[i]]<-W_test_t_Fixed_10[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_10<-as.list(p.adjust(p_values_Fixed_10,method="bonferroni"))
B_Significative_Fixed_10<-list()
for (i in 1:20){B_Significative_Fixed_10[[i]]<-
  if (B_corrected_p_values_Fixed_10[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_10<-as.list(p.adjust(p_values_Fixed_10,method="BH"))
BH_Significative_Fixed_10<-list()
for (i in 1:20){BH_Significative_Fixed_10[[i]]<-
  if (BH_corrected_p_values_Fixed_10[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_10<-paste(p_values_Fixed_10, collapse = NULL)
Statistic_W_Fixed_10<-paste(Statistic_W_Fixed_10, collapse = NULL)
B_corrected_p_values_Fixed_10<-paste(B_corrected_p_values_Fixed_10, collapse = NULL)
B_Significative_Fixed_10<-paste(B_Significative_Fixed_10, collapse = NULL)
BH_corrected_p_values_Fixed_10<-paste(BH_corrected_p_values_Fixed_10, collapse = NULL)
BH_Significative_Fixed_10<-paste(BH_Significative_Fixed_10, collapse = NULL)

#Results
Length_groups<-c("10 vs 1","10 vs 2","10 vs 3","10 vs 4","10 vs 5","10 vs 6","10 vs 7","10 vs 8","10 vs 9","10 vs 10","10 vs 11","10 vs 12","10 vs 13","10 vs 14","10 vs 15","10 vs 16","10 vs 17","10 vs 18","10 vs 19","10 vs 20")
Fixed_group10_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_10,Statistic_W_Fixed_10,B_corrected_p_values_Fixed_10,
                                                          B_Significative_Fixed_10,BH_corrected_p_values_Fixed_10,BH_Significative_Fixed_10)
write.xlsx(Fixed_group10_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_10_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------11----------------------------------------------------------------------------
W_test_t_Fixed_11<-list()
for (i in 1:20){W_test_t_Fixed_11[[i]]<-wilcox.test(Fixed_exp_list[[11]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_11<-list()
for (i in 1:20){p_values_Fixed_11[[i]]<-W_test_t_Fixed_11[[i]][["p.value"]]}
Statistic_W_Fixed_11<-list()
for (i in 1:20){Statistic_W_Fixed_11[[i]]<-W_test_t_Fixed_11[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_11<-as.list(p.adjust(p_values_Fixed_11,method="bonferroni"))
B_Significative_Fixed_11<-list()
for (i in 1:20){B_Significative_Fixed_11[[i]]<-
  if (B_corrected_p_values_Fixed_11[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_11<-as.list(p.adjust(p_values_Fixed_11,method="BH"))
BH_Significative_Fixed_11<-list()
for (i in 1:20){BH_Significative_Fixed_11[[i]]<-
  if (BH_corrected_p_values_Fixed_11[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_11<-paste(p_values_Fixed_11, collapse = NULL)
Statistic_W_Fixed_11<-paste(Statistic_W_Fixed_11, collapse = NULL)
B_corrected_p_values_Fixed_11<-paste(B_corrected_p_values_Fixed_11, collapse = NULL)
B_Significative_Fixed_11<-paste(B_Significative_Fixed_11, collapse = NULL)
BH_corrected_p_values_Fixed_11<-paste(BH_corrected_p_values_Fixed_11, collapse = NULL)
BH_Significative_Fixed_11<-paste(BH_Significative_Fixed_11, collapse = NULL)

#Results
Length_groups<-c("11 vs 1","11 vs 2","11 vs 3","11 vs 4","11 vs 5","11 vs 6","11 vs 7","11 vs 8","11 vs 9","11 vs 10","11 vs 11","11 vs 12","11 vs 13","11 vs 14","11 vs 15","11 vs 16","11 vs 17","11 vs 18","11 vs 19","11 vs 20")
Fixed_group11_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_11,Statistic_W_Fixed_11,B_corrected_p_values_Fixed_11,
                                                          B_Significative_Fixed_11,BH_corrected_p_values_Fixed_11,BH_Significative_Fixed_11)
write.xlsx(Fixed_group11_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_11_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------12----------------------------------------------------------------------------
W_test_t_Fixed_12<-list()
for (i in 1:20){W_test_t_Fixed_12[[i]]<-wilcox.test(Fixed_exp_list[[12]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_12<-list()
for (i in 1:20){p_values_Fixed_12[[i]]<-W_test_t_Fixed_12[[i]][["p.value"]]}
Statistic_W_Fixed_12<-list()
for (i in 1:20){Statistic_W_Fixed_12[[i]]<-W_test_t_Fixed_12[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_12<-as.list(p.adjust(p_values_Fixed_12,method="bonferroni"))
B_Significative_Fixed_12<-list()
for (i in 1:20){B_Significative_Fixed_12[[i]]<-
  if (B_corrected_p_values_Fixed_12[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_12<-as.list(p.adjust(p_values_Fixed_12,method="BH"))
BH_Significative_Fixed_12<-list()
for (i in 1:20){BH_Significative_Fixed_12[[i]]<-
  if (BH_corrected_p_values_Fixed_12[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_12<-paste(p_values_Fixed_12, collapse = NULL)
Statistic_W_Fixed_12<-paste(Statistic_W_Fixed_12, collapse = NULL)
B_corrected_p_values_Fixed_12<-paste(B_corrected_p_values_Fixed_12, collapse = NULL)
B_Significative_Fixed_12<-paste(B_Significative_Fixed_12, collapse = NULL)
BH_corrected_p_values_Fixed_12<-paste(BH_corrected_p_values_Fixed_12, collapse = NULL)
BH_Significative_Fixed_12<-paste(BH_Significative_Fixed_12, collapse = NULL)

#Results
Length_groups<-c("12 vs 1","12 vs 2","12 vs 3","12 vs 4","12 vs 5","12 vs 6","12 vs 7","12 vs 8","12 vs 9","12 vs 10","12 vs 11","12 vs 12","12 vs 13","12 vs 14","12 vs 15","12 vs 16","12 vs 17","12 vs 18","12 vs 19","12 vs 20")
Fixed_group12_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_12,Statistic_W_Fixed_12,B_corrected_p_values_Fixed_12,
                                                          B_Significative_Fixed_12,BH_corrected_p_values_Fixed_12,BH_Significative_Fixed_12)
write.xlsx(Fixed_group12_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_12_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)
#-----------------------------------------------------------------13----------------------------------------------------------------------------
W_test_t_Fixed_13<-list()
for (i in 1:20){W_test_t_Fixed_13[[i]]<-wilcox.test(Fixed_exp_list[[13]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_13<-list()
for (i in 1:20){p_values_Fixed_13[[i]]<-W_test_t_Fixed_13[[i]][["p.value"]]}
Statistic_W_Fixed_13<-list()
for (i in 1:20){Statistic_W_Fixed_13[[i]]<-W_test_t_Fixed_13[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_13<-as.list(p.adjust(p_values_Fixed_13,method="bonferroni"))
B_Significative_Fixed_13<-list()
for (i in 1:20){B_Significative_Fixed_13[[i]]<-
  if (B_corrected_p_values_Fixed_13[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_13<-as.list(p.adjust(p_values_Fixed_13,method="BH"))
BH_Significative_Fixed_13<-list()
for (i in 1:20){BH_Significative_Fixed_13[[i]]<-
  if (BH_corrected_p_values_Fixed_13[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_13<-paste(p_values_Fixed_13, collapse = NULL)
Statistic_W_Fixed_13<-paste(Statistic_W_Fixed_13, collapse = NULL)
B_corrected_p_values_Fixed_13<-paste(B_corrected_p_values_Fixed_13, collapse = NULL)
B_Significative_Fixed_13<-paste(B_Significative_Fixed_13, collapse = NULL)
BH_corrected_p_values_Fixed_13<-paste(BH_corrected_p_values_Fixed_13, collapse = NULL)
BH_Significative_Fixed_13<-paste(BH_Significative_Fixed_13, collapse = NULL)

#Results
Length_groups<-c("13 vs 1","13 vs 2","13 vs 3","13 vs 4","13 vs 5","13 vs 6","13 vs 7","13 vs 8","13 vs 9","13 vs 10","13 vs 11","13 vs 12","13 vs 13","13 vs 14","13 vs 15","13 vs 16","13 vs 17","13 vs 18","13 vs 19","13 vs 20")
Fixed_group13_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_13,Statistic_W_Fixed_13,B_corrected_p_values_Fixed_13,
                                                          B_Significative_Fixed_13,BH_corrected_p_values_Fixed_13,BH_Significative_Fixed_13)
write.xlsx(Fixed_group13_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_13_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)
#-----------------------------------------------------------------14----------------------------------------------------------------------------
W_test_t_Fixed_14<-list()
for (i in 1:20){W_test_t_Fixed_14[[i]]<-wilcox.test(Fixed_exp_list[[14]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_14<-list()
for (i in 1:20){p_values_Fixed_14[[i]]<-W_test_t_Fixed_14[[i]][["p.value"]]}
Statistic_W_Fixed_14<-list()
for (i in 1:20){Statistic_W_Fixed_14[[i]]<-W_test_t_Fixed_14[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_14<-as.list(p.adjust(p_values_Fixed_14,method="bonferroni"))
B_Significative_Fixed_14<-list()
for (i in 1:20){B_Significative_Fixed_14[[i]]<-
  if (B_corrected_p_values_Fixed_14[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_14<-as.list(p.adjust(p_values_Fixed_14,method="BH"))
BH_Significative_Fixed_14<-list()
for (i in 1:20){BH_Significative_Fixed_14[[i]]<-
  if (BH_corrected_p_values_Fixed_14[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_14<-paste(p_values_Fixed_14, collapse = NULL)
Statistic_W_Fixed_14<-paste(Statistic_W_Fixed_14, collapse = NULL)
B_corrected_p_values_Fixed_14<-paste(B_corrected_p_values_Fixed_14, collapse = NULL)
B_Significative_Fixed_14<-paste(B_Significative_Fixed_14, collapse = NULL)
BH_corrected_p_values_Fixed_14<-paste(BH_corrected_p_values_Fixed_14, collapse = NULL)
BH_Significative_Fixed_14<-paste(BH_Significative_Fixed_14, collapse = NULL)

#Results
Length_groups<-c("14 vs 1","14 vs 2","14 vs 3","14 vs 4","14 vs 5","14 vs 6","14 vs 7","14 vs 8","14 vs 9","14 vs 10","14 vs 11","14 vs 12","14 vs 13","14 vs 14","14 vs 15","14 vs 16","14 vs 17","14 vs 18","14 vs 19","14 vs 20")
Fixed_group14_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_14,Statistic_W_Fixed_14,B_corrected_p_values_Fixed_14,
                                                          B_Significative_Fixed_14,BH_corrected_p_values_Fixed_14,BH_Significative_Fixed_14)
write.xlsx(Fixed_group14_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_14_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------15----------------------------------------------------------------------------
W_test_t_Fixed_15<-list()
for (i in 1:20){W_test_t_Fixed_15[[i]]<-wilcox.test(Fixed_exp_list[[15]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_15<-list()
for (i in 1:20){p_values_Fixed_15[[i]]<-W_test_t_Fixed_15[[i]][["p.value"]]}
Statistic_W_Fixed_15<-list()
for (i in 1:20){Statistic_W_Fixed_15[[i]]<-W_test_t_Fixed_15[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_15<-as.list(p.adjust(p_values_Fixed_15,method="bonferroni"))
B_Significative_Fixed_15<-list()
for (i in 1:20){B_Significative_Fixed_15[[i]]<-
  if (B_corrected_p_values_Fixed_15[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_15<-as.list(p.adjust(p_values_Fixed_15,method="BH"))
BH_Significative_Fixed_15<-list()
for (i in 1:20){BH_Significative_Fixed_15[[i]]<-
  if (BH_corrected_p_values_Fixed_15[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_15<-paste(p_values_Fixed_15, collapse = NULL)
Statistic_W_Fixed_15<-paste(Statistic_W_Fixed_15, collapse = NULL)
B_corrected_p_values_Fixed_15<-paste(B_corrected_p_values_Fixed_15, collapse = NULL)
B_Significative_Fixed_15<-paste(B_Significative_Fixed_15, collapse = NULL)
BH_corrected_p_values_Fixed_15<-paste(BH_corrected_p_values_Fixed_15, collapse = NULL)
BH_Significative_Fixed_15<-paste(BH_Significative_Fixed_15, collapse = NULL)

#Results
Length_groups<-c("15 vs 1","15 vs 2","15 vs 3","15 vs 4","15 vs 5","15 vs 6","15 vs 7","15 vs 8","15 vs 9","15 vs 10","15 vs 11","15 vs 12","15 vs 13","15 vs 14","15 vs 15","15 vs 16","15 vs 17","15 vs 18","15 vs 19","15 vs 20")
Fixed_group15_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_15,Statistic_W_Fixed_15,B_corrected_p_values_Fixed_15,
                                                          B_Significative_Fixed_15,BH_corrected_p_values_Fixed_15,BH_Significative_Fixed_15)
write.xlsx(Fixed_group15_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_15_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)
#-----------------------------------------------------------------16----------------------------------------------------------------------------
W_test_t_Fixed_16<-list()
for (i in 1:20){W_test_t_Fixed_16[[i]]<-wilcox.test(Fixed_exp_list[[16]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_16<-list()
for (i in 1:20){p_values_Fixed_16[[i]]<-W_test_t_Fixed_16[[i]][["p.value"]]}
Statistic_W_Fixed_16<-list()
for (i in 1:20){Statistic_W_Fixed_16[[i]]<-W_test_t_Fixed_16[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_16<-as.list(p.adjust(p_values_Fixed_16,method="bonferroni"))
B_Significative_Fixed_16<-list()
for (i in 1:20){B_Significative_Fixed_16[[i]]<-
  if (B_corrected_p_values_Fixed_16[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_16<-as.list(p.adjust(p_values_Fixed_16,method="BH"))
BH_Significative_Fixed_16<-list()
for (i in 1:20){BH_Significative_Fixed_16[[i]]<-
  if (BH_corrected_p_values_Fixed_16[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_16<-paste(p_values_Fixed_16, collapse = NULL)
Statistic_W_Fixed_16<-paste(Statistic_W_Fixed_16, collapse = NULL)
B_corrected_p_values_Fixed_16<-paste(B_corrected_p_values_Fixed_16, collapse = NULL)
B_Significative_Fixed_16<-paste(B_Significative_Fixed_16, collapse = NULL)
BH_corrected_p_values_Fixed_16<-paste(BH_corrected_p_values_Fixed_16, collapse = NULL)
BH_Significative_Fixed_16<-paste(BH_Significative_Fixed_16, collapse = NULL)

#Results
Length_groups<-c("16 vs 1","16 vs 2","16 vs 3","16 vs 4","16 vs 5","16 vs 6","16 vs 7","16 vs 8","16 vs 9","16 vs 10","16 vs 11","16 vs 12","16 vs 13","16 vs 14","16 vs 15","16 vs 16","16 vs 17","16 vs 18","16 vs 19","16 vs 20")
Fixed_group16_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_16,Statistic_W_Fixed_16,B_corrected_p_values_Fixed_16,
                                                          B_Significative_Fixed_16,BH_corrected_p_values_Fixed_16,BH_Significative_Fixed_16)
write.xlsx(Fixed_group16_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_16_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------17----------------------------------------------------------------------------
W_test_t_Fixed_17<-list()
for (i in 1:20){W_test_t_Fixed_17[[i]]<-wilcox.test(Fixed_exp_list[[17]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_17<-list()
for (i in 1:20){p_values_Fixed_17[[i]]<-W_test_t_Fixed_17[[i]][["p.value"]]}
Statistic_W_Fixed_17<-list()
for (i in 1:20){Statistic_W_Fixed_17[[i]]<-W_test_t_Fixed_17[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_17<-as.list(p.adjust(p_values_Fixed_17,method="bonferroni"))
B_Significative_Fixed_17<-list()
for (i in 1:20){B_Significative_Fixed_17[[i]]<-
  if (B_corrected_p_values_Fixed_17[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_17<-as.list(p.adjust(p_values_Fixed_17,method="BH"))
BH_Significative_Fixed_17<-list()
for (i in 1:20){BH_Significative_Fixed_17[[i]]<-
  if (BH_corrected_p_values_Fixed_17[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_17<-paste(p_values_Fixed_17, collapse = NULL)
Statistic_W_Fixed_17<-paste(Statistic_W_Fixed_17, collapse = NULL)
B_corrected_p_values_Fixed_17<-paste(B_corrected_p_values_Fixed_17, collapse = NULL)
B_Significative_Fixed_17<-paste(B_Significative_Fixed_17, collapse = NULL)
BH_corrected_p_values_Fixed_17<-paste(BH_corrected_p_values_Fixed_17, collapse = NULL)
BH_Significative_Fixed_17<-paste(BH_Significative_Fixed_17, collapse = NULL)

#Results
Length_groups<-c("17 vs 1","17 vs 2","17 vs 3","17 vs 4","17 vs 5","17 vs 6","17 vs 7","17 vs 8","17 vs 9","17 vs 10","17 vs 11","17 vs 12","17 vs 13","17 vs 14","17 vs 15","17 vs 16","17 vs 17","17 vs 18","17 vs 19","17 vs 20")
Fixed_group17_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_17,Statistic_W_Fixed_17,B_corrected_p_values_Fixed_17,
                                                          B_Significative_Fixed_17,BH_corrected_p_values_Fixed_17,BH_Significative_Fixed_17)
write.xlsx(Fixed_group17_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_17_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)


#-----------------------------------------------------------------18----------------------------------------------------------------------------
W_test_t_Fixed_18<-list()
for (i in 1:20){W_test_t_Fixed_18[[i]]<-wilcox.test(Fixed_exp_list[[18]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_18<-list()
for (i in 1:20){p_values_Fixed_18[[i]]<-W_test_t_Fixed_18[[i]][["p.value"]]}
Statistic_W_Fixed_18<-list()
for (i in 1:20){Statistic_W_Fixed_18[[i]]<-W_test_t_Fixed_18[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_18<-as.list(p.adjust(p_values_Fixed_18,method="bonferroni"))
B_Significative_Fixed_18<-list()
for (i in 1:20){B_Significative_Fixed_18[[i]]<-
  if (B_corrected_p_values_Fixed_18[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_18<-as.list(p.adjust(p_values_Fixed_18,method="BH"))
BH_Significative_Fixed_18<-list()
for (i in 1:20){BH_Significative_Fixed_18[[i]]<-
  if (BH_corrected_p_values_Fixed_18[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_18<-paste(p_values_Fixed_18, collapse = NULL)
Statistic_W_Fixed_18<-paste(Statistic_W_Fixed_18, collapse = NULL)
B_corrected_p_values_Fixed_18<-paste(B_corrected_p_values_Fixed_18, collapse = NULL)
B_Significative_Fixed_18<-paste(B_Significative_Fixed_18, collapse = NULL)
BH_corrected_p_values_Fixed_18<-paste(BH_corrected_p_values_Fixed_18, collapse = NULL)
BH_Significative_Fixed_18<-paste(BH_Significative_Fixed_18, collapse = NULL)

#Results
Length_groups<-c("18 vs 1","18 vs 2","18 vs 3","18 vs 4","18 vs 5","18 vs 6","18 vs 7","18 vs 8","18 vs 9","18 vs 10","18 vs 11","18 vs 12","18 vs 13","18 vs 14","18 vs 15","18 vs 16","18 vs 17","18 vs 18","18 vs 19","18 vs 20")
Fixed_group18_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_18,Statistic_W_Fixed_18,B_corrected_p_values_Fixed_18,
                                                          B_Significative_Fixed_18,BH_corrected_p_values_Fixed_18,BH_Significative_Fixed_18)
write.xlsx(Fixed_group18_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_18_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)

#-----------------------------------------------------------------19----------------------------------------------------------------------------
W_test_t_Fixed_19<-list()
for (i in 1:20){W_test_t_Fixed_19[[i]]<-wilcox.test(Fixed_exp_list[[19]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_19<-list()
for (i in 1:20){p_values_Fixed_19[[i]]<-W_test_t_Fixed_19[[i]][["p.value"]]}
Statistic_W_Fixed_19<-list()
for (i in 1:20){Statistic_W_Fixed_19[[i]]<-W_test_t_Fixed_19[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_19<-as.list(p.adjust(p_values_Fixed_19,method="bonferroni"))
B_Significative_Fixed_19<-list()
for (i in 1:20){B_Significative_Fixed_19[[i]]<-
  if (B_corrected_p_values_Fixed_19[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_19<-as.list(p.adjust(p_values_Fixed_19,method="BH"))
BH_Significative_Fixed_19<-list()
for (i in 1:20){BH_Significative_Fixed_19[[i]]<-
  if (BH_corrected_p_values_Fixed_19[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_19<-paste(p_values_Fixed_19, collapse = NULL)
Statistic_W_Fixed_19<-paste(Statistic_W_Fixed_19, collapse = NULL)
B_corrected_p_values_Fixed_19<-paste(B_corrected_p_values_Fixed_19, collapse = NULL)
B_Significative_Fixed_19<-paste(B_Significative_Fixed_19, collapse = NULL)
BH_corrected_p_values_Fixed_19<-paste(BH_corrected_p_values_Fixed_19, collapse = NULL)
BH_Significative_Fixed_19<-paste(BH_Significative_Fixed_19, collapse = NULL)

#Results
Length_groups<-c("19 vs 1","19 vs 2","19 vs 3","19 vs 4","19 vs 5","19 vs 6","19 vs 7","19 vs 8","19 vs 9","19 vs 10","19 vs 11","19 vs 12","19 vs 13","19 vs 14","19 vs 15","19 vs 16","19 vs 17","19 vs 18","19 vs 19","19 vs 20")
Fixed_group19_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_19,Statistic_W_Fixed_19,B_corrected_p_values_Fixed_19,
                                                          B_Significative_Fixed_19,BH_corrected_p_values_Fixed_19,BH_Significative_Fixed_19)
write.xlsx(Fixed_group19_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_19_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)
#-----------------------------------------------------------------20----------------------------------------------------------------------------
W_test_t_Fixed_20<-list()
for (i in 1:20){W_test_t_Fixed_20[[i]]<-wilcox.test(Fixed_exp_list[[20]],Fixed_exp_list[[i]], paired=FALSE)}
p_values_Fixed_20<-list()
for (i in 1:20){p_values_Fixed_20[[i]]<-W_test_t_Fixed_20[[i]][["p.value"]]}
Statistic_W_Fixed_20<-list()
for (i in 1:20){Statistic_W_Fixed_20[[i]]<-W_test_t_Fixed_20[[i]][["statistic"]]}


#p-values adjustment
B_corrected_p_values_Fixed_20<-as.list(p.adjust(p_values_Fixed_20,method="bonferroni"))
B_Significative_Fixed_20<-list()
for (i in 1:20){B_Significative_Fixed_20[[i]]<-
  if (B_corrected_p_values_Fixed_20[[i]]<=0.05){"YES"} else {"NO"}}
BH_corrected_p_values_Fixed_20<-as.list(p.adjust(p_values_Fixed_20,method="BH"))
BH_Significative_Fixed_20<-list()
for (i in 1:20){BH_Significative_Fixed_20[[i]]<-
  if (BH_corrected_p_values_Fixed_20[[i]]<=0.05){"YES"} else {"NO"}}

p_values_Fixed_20<-paste(p_values_Fixed_20, collapse = NULL)
Statistic_W_Fixed_20<-paste(Statistic_W_Fixed_20, collapse = NULL)
B_corrected_p_values_Fixed_20<-paste(B_corrected_p_values_Fixed_20, collapse = NULL)
B_Significative_Fixed_20<-paste(B_Significative_Fixed_20, collapse = NULL)
BH_corrected_p_values_Fixed_20<-paste(BH_corrected_p_values_Fixed_20, collapse = NULL)
BH_Significative_Fixed_20<-paste(BH_Significative_Fixed_20, collapse = NULL)

#Results
Length_groups<-c("20 vs 1","20 vs 2","20 vs 3","20 vs 4","20 vs 5","20 vs 6","20 vs 7","20 vs 8","20 vs 9","20 vs 10","20 vs 11","20 vs 12","20 vs 13","20 vs 14","20 vs 15","20 vs 16","20 vs 17","20 vs 18","20 vs 19","20 vs 20")
Fixed_group20_against_rest_statistics_summary<-data.frame(Length_groups,p_values_Fixed_20,Statistic_W_Fixed_20,B_corrected_p_values_Fixed_20,
                                                          B_Significative_Fixed_20,BH_corrected_p_values_Fixed_20,BH_Significative_Fixed_20)
write.xlsx(Fixed_group20_against_rest_statistics_summary, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/",paste(time.name),"Length_groups_statistics_summary.xlsx"))
           , sheetName = "Fixed_20_vs_rest", col.names = TRUE, row.names = TRUE, append = TRUE)


#Plots to represent statistics
#Load
Sorted_Effect_methanol_depending_length_to_plot<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Sorted_Effect_methanol_depending_length_to_plot.rds")
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
saveRDS(Sorted_Effect_methanol_depending_length_to_plot, file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Sorted_Effect_methanol_depending_length_to_plot_2.rds")
#Load
Sorted_Effect_methanol_depending_length_to_plot<- readRDS(file="F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/11 Evaluate effect depending on trancript length/R objects/Sorted_Effect_methanol_depending_length_to_plot_2.rds")

ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=reorder(Groups_range,Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Average expression level as a function of transcript length")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+ 
  ylim(0,2)+
  xlab("Transcript length categories") + ylab("Average expression")

ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=reorder(Groups_range,Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA,data=subset(Sorted_Effect_methanol_depending_length_to_plot,Condition=="Fresh"))+labs(title="Average expression level as a function of transcript length")+
  scale_fill_manual(values=c("#F8766D"))+ 
  ylim(0,2)+
  xlab("Transcript length categories") + ylab("Average expression")

ggplot(Sorted_Effect_methanol_depending_length_to_plot, aes(x=reorder(Groups_range,Groups), y=Levels_of_expression, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA,data=subset(Sorted_Effect_methanol_depending_length_to_plot,Condition=="Fixed"))+labs(title="Average expression level as a function of transcript length")+
  scale_fill_manual(values=c("#00BFC4"))+ 
  ylim(0,2)+
  xlab("Transcript length categories") + ylab("Average expression")

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

ggplot(Ratio_df, aes(x=reorder(Groups_range,Groups), y=Ratio)) + 
  geom_boxplot(outlier.shape = NA)+labs(title="Ratio (1+Fresh/1+Fixed) Average expression level as a function of transcript length")+
  scale_fill_manual(values=c("#9F9F9F"))+ 
  ylim(0.625,1.75)+
  xlab("Transcript length categories") + ylab("Ratio 1+average expression Fresh/Fixed")+ theme_classic()


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



