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
install.packages('PupillometryR')
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
#library(ggtree)
library(patchwork) 
library(PupillometryR)
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


#chack the mitochondrial genes that we have in our data set
Row_names<-row.names(pilot.data)
Row_names
Mito_genes_names<-grep("^mt",Row_names, value=TRUE)
Mito_genes_names

Mito_genes_names_2<-grep("^mt-",Row_names, value=TRUE)
Mito_genes_names_2


#when we create the seurat object we need to ask to remove the transcripts belonging to the ercc matrix and 
#add total_ercc_reads and percent.ERCC to object@meta.data in the total_ercc_reads and percent.ERCC column respectively
pilot <- CreateSeuratObject(counts = pilot.data[-ercc,], project = "Pilot analysis", min.cells = 0, min.features = 0)
pilot


pilot$orig.ident <- factor(pilot$orig.ident, levels=c("JEN-AU-s004","JEN-AU-s002"))

current.ids <- c("JEN-AU-s004","JEN-AU-s002")
new.ids <- c("Fresh_cells","Fixed_cells")
pilot@active.ident <- plyr::mapvalues(x = pilot@active.ident, from = current.ids, to = new.ids)
pilot$Condition <- Idents(pilot)

pilot$Condition <- factor(x=pilot$Condition, levels=c("Fresh_cells","Fixed_cells"))

table(pilot$Condition)
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
pilot[["percent_cyto"]] <- 100-(pilot$percent_mt)
pilot[["percent_ERCC"]] <- PercentageFeatureSet(pilot, pattern = "^ERCC")
pilot[["log10_nCount_RNA"]] <- log10(pilot$nCount_RNA)
pilot[["log10_nFeature_RNA"]] <- log10(pilot$nFeature_RNA)
pilot[["percent_ERCC_after_removal"]] <- PercentageFeatureSet(pilot, pattern = "^ERCC")
pilot[["ratio_mol_gen"]]<-pilot$nCount_RNA/pilot$nFeature_RNA
Control_of_ERCC_removal<- summary(pilot$percent_ERCC_after_removal)
Control_of_ERCC_removal
Control_percent_mt<-summary(pilot$percent_mt)
Control_percent_mt



pilot.list<-SplitObject(pilot, split.by = "Condition")

Fresh_cells<-pilot.list$Fresh_cells
Fixed_cells<-pilot.list$Fixed_cells


#Density plots
log10_nCount_RNA<-pilot$log10_nCount_RNA
log10_nFeature_RNA<-pilot$log10_nFeature_RNA
Condition<-pilot$Condition
percent_mt<-pilot$percent_mt
Density_plot_df<-data.frame(Condition,log10_nCount_RNA,log10_nFeature_RNA,percent_mt)

ggplot(Density_plot_df, aes(x=log10_nCount_RNA, color=Condition)) +
  geom_density(size = 1.5)+scale_color_manual(values=c("#F8766D", "#00BFC4")) +
  geom_vline(xintercept =2.92942) + geom_vline(xintercept =4.54407) + theme_classic()

ggplot(Density_plot_df, aes(x=log10_nFeature_RNA, color=Condition)) +
  geom_density(size = 1.5, adjust = 0.05)+scale_color_manual(values=c("#F8766D", "#00BFC4")) +
  geom_vline(xintercept =2.69897) + theme_classic()
  #+xlim(2.7,3.9)

ggplot(Density_plot_df, aes(x=percent_mt, color=Condition)) +
  geom_density(size = 1.5)+scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()

#Density plots with box plots_log10
#Arrangement of the dataframe
Density_plot_df$Cells<-rownames(Density_plot_df)
Density__and_box_plot_df_log10_nCount_RNA<-Density_plot_df[,-c(3,4)]
Density__and_box_plot_df_log10_nCount_RNA$QC_Parameter<-rep("log10_nCount_RNA", times=768)
ggplot(Density__and_box_plot_df_log10_nCount_RNA, aes(x=QC_Parameter, y=log10_nCount_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =2.92942) + geom_hline(yintercept =4.54407)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)
  

Density__and_box_plot_df_log10_nFeature_RNA<-Density_plot_df[,-c(2,4)]
Density__and_box_plot_df_log10_nFeature_RNA$QC_Parameter<-rep("log10_nFeature_RNA", times=768)
ggplot(Density__and_box_plot_df_log10_nFeature_RNA, aes(x=QC_Parameter, y=log10_nFeature_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =2.69897)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)

Density__and_box_plot_df_percent_mt<-Density_plot_df[,-c(3,2)]
Density__and_box_plot_df_percent_mt$QC_Parameter<-rep("percent_mt", times=768)
ggplot(Density__and_box_plot_df_percent_mt, aes(x=QC_Parameter, y=percent_mt, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)



#Density plots with box plots_NO log10!!!!
#Arrangement of the dataframe
nCount_RNA<-pilot$nCount_RNA
nFeature_RNA<-pilot$nFeature_RNA
Condition<-pilot$Condition
percent_mt<-pilot$percent_mt
Density_plot_df<-data.frame(Condition,nCount_RNA,nFeature_RNA,percent_mt)
Density_plot_df$Cells<-rownames(Density_plot_df)

Density__and_box_plot_df_nCount_RNA<-Density_plot_df[,-c(3,4)]
Density__and_box_plot_df_nCount_RNA$QC_Parameter<-rep("nCount_RNA", times=768)
ggplot(Density__and_box_plot_df_nCount_RNA, aes(x=QC_Parameter, y=nCount_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =800) + geom_hline(yintercept =35000)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)


Density__and_box_plot_df_nFeature_RNA<-Density_plot_df[,-c(2,4)]
Density__and_box_plot_df_nFeature_RNA$QC_Parameter<-rep("nFeature_RNA", times=768)
ggplot(Density__and_box_plot_df_nFeature_RNA, aes(x=QC_Parameter, y=nFeature_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =500)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)

Density__and_box_plot_df_percent_mt<-Density_plot_df[,-c(3,2)]
Density__and_box_plot_df_percent_mt$QC_Parameter<-rep("percent_mt", times=768)
ggplot(Density__and_box_plot_df_percent_mt, aes(x=QC_Parameter, y=percent_mt, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F) + ylim(0, 100)











#--------------------------------------------------------filtering---------------------------------------------------------------------
Fresh_cells_filt<-subset(Fresh_cells, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
Fixed_cells_filt<-subset(Fixed_cells, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
pilot_filt<-subset(pilot, subset = nCount_RNA > 800 & nCount_RNA < 35000 & percent_mt < 100 & 
                      nFeature_RNA >500 & total_ercc_reads >500 & ratio_mol_gen>1.2)
Fresh_cells_filt
Fixed_cells_filt
pilot_filt

#Density plots
log10_nCount_RNA<-pilot_filt$log10_nCount_RNA
log10_nFeature_RNA<-pilot_filt$log10_nFeature_RNA
Condition<-pilot_filt$Condition
percent_mt<-pilot_filt$percent_mt
Density_plot_df<-data.frame(Condition,log10_nCount_RNA,log10_nFeature_RNA,percent_mt)

ggplot(Density_plot_df, aes(x=log10_nCount_RNA, color=Condition)) +
  geom_density(size = 1.5)+scale_color_manual(values=c("#F8766D", "#00BFC4")) +
  geom_vline(xintercept =2.92942) + geom_vline(xintercept =4.54407) + theme_classic()

ggplot(Density_plot_df, aes(x=log10_nFeature_RNA, color=Condition)) +
  geom_density(size = 1.5)+scale_color_manual(values=c("#F8766D", "#00BFC4")) +
  geom_vline(xintercept =2.69897) + theme_classic()

ggplot(Density_plot_df, aes(x=percent_mt, color=Condition)) +
  geom_density(size = 1.5)+scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()



#Density plots with box plots
#Arrangement of the dataframe
nCount_RNA<-pilot_filt$nCount_RNA
nFeature_RNA<-pilot_filt$nFeature_RNA
Condition<-pilot_filt$Condition
percent_mt<-pilot_filt$percent_mt
Density_plot_df<-data.frame(Condition,nCount_RNA,nFeature_RNA,percent_mt)
Density_plot_df$Cells<-rownames(Density_plot_df)

Density__and_box_plot_df_nCount_RNA<-Density_plot_df[,-c(3,4)]
Density__and_box_plot_df_nCount_RNA$QC_Parameter<-rep("nCount_RNA", times=196)
ggplot(Density__and_box_plot_df_nCount_RNA, aes(x=QC_Parameter, y=nCount_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =800) + geom_hline(yintercept =35000)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)


Density__and_box_plot_df_nFeature_RNA<-Density_plot_df[,-c(2,4)]
Density__and_box_plot_df_nFeature_RNA$QC_Parameter<-rep("nFeature_RNA", times=196)
ggplot(Density__and_box_plot_df_nFeature_RNA, aes(x=QC_Parameter, y=nFeature_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =500)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)

Density__and_box_plot_df_percent_mt<-Density_plot_df[,-c(3,2)]
Density__and_box_plot_df_percent_mt$QC_Parameter<-rep("percent_mt", times=196)
ggplot(Density__and_box_plot_df_percent_mt, aes(x=QC_Parameter, y=percent_mt, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F) + ylim(0, 100)


#Density plots with box plots_log10
#Arrangement of the dataframe
#Density plots
log10_nCount_RNA<-pilot_filt$log10_nCount_RNA
log10_nFeature_RNA<-pilot_filt$log10_nFeature_RNA
Condition<-pilot_filt$Condition
percent_mt<-pilot_filt$percent_mt
Density_plot_df<-data.frame(Condition,log10_nCount_RNA,log10_nFeature_RNA,percent_mt)

Density_plot_df$Cells<-rownames(Density_plot_df)
Density__and_box_plot_df_log10_nCount_RNA<-Density_plot_df[,-c(3,4)]
Density__and_box_plot_df_log10_nCount_RNA$QC_Parameter<-rep("log10_nCount_RNA", times=196)
ggplot(Density__and_box_plot_df_log10_nCount_RNA, aes(x=QC_Parameter, y=log10_nCount_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =2.92942) + geom_hline(yintercept =4.54407)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)


Density__and_box_plot_df_log10_nFeature_RNA<-Density_plot_df[,-c(2,4)]
Density__and_box_plot_df_log10_nFeature_RNA$QC_Parameter<-rep("log10_nFeature_RNA", times=196)
ggplot(Density__and_box_plot_df_log10_nFeature_RNA, aes(x=QC_Parameter, y=log10_nFeature_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =2.69897)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)

Density__and_box_plot_df_percent_mt<-Density_plot_df[,-c(3,2)]
Density__and_box_plot_df_percent_mt$QC_Parameter<-rep("percent_mt", times=196)
ggplot(Density__and_box_plot_df_percent_mt, aes(x=QC_Parameter, y=percent_mt, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)+ylim(0,100)



#--------------------------------------------------------filtering with only n_Features threshold---------------------------------------------------------------------
Fresh_cells_filt<-subset(Fresh_cells, subset = nFeature_RNA >500)
Fixed_cells_filt<-subset(Fixed_cells, subset = nFeature_RNA >500)
pilot_filt<-subset(pilot, subset = nFeature_RNA >500)
Fresh_cells_filt
Fixed_cells_filt
pilot_filt


#Density plots
log10_nCount_RNA<-pilot_filt$log10_nCount_RNA
log10_nFeature_RNA<-pilot_filt$log10_nFeature_RNA
Condition<-pilot_filt$Condition
percent_mt<-pilot_filt$percent_mt
Density_plot_df<-data.frame(Condition,log10_nCount_RNA,log10_nFeature_RNA,percent_mt)

ggplot(Density_plot_df, aes(x=log10_nCount_RNA, color=Condition)) +
  geom_density(size = 1.5)+scale_color_manual(values=c("#F8766D", "#00BFC4")) +
  geom_vline(xintercept =2.92942) + geom_vline(xintercept =4.54407) + theme_classic()

ggplot(Density_plot_df, aes(x=log10_nFeature_RNA, color=Condition)) +
  geom_density(size = 1.5)+scale_color_manual(values=c("#F8766D", "#00BFC4")) +
  geom_vline(xintercept =2.69897) + theme_classic()

ggplot(Density_plot_df, aes(x=percent_mt, color=Condition)) +
  geom_density(size = 1.5)+scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()




#Density plots with box plots
#Arrangement of the dataframe
nCount_RNA<-pilot_filt$nCount_RNA
nFeature_RNA<-pilot_filt$nFeature_RNA
Condition<-pilot_filt$Condition
percent_mt<-pilot_filt$percent_mt
Density_plot_df<-data.frame(Condition,nCount_RNA,nFeature_RNA,percent_mt)
Density_plot_df$Cells<-rownames(Density_plot_df)

Density__and_box_plot_df_nCount_RNA<-Density_plot_df[,-c(3,4)]
Density__and_box_plot_df_nCount_RNA$QC_Parameter<-rep("nCount_RNA", times=197)
ggplot(Density__and_box_plot_df_nCount_RNA, aes(x=QC_Parameter, y=nCount_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =800) + geom_hline(yintercept =35000)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)


Density__and_box_plot_df_nFeature_RNA<-Density_plot_df[,-c(2,4)]
Density__and_box_plot_df_nFeature_RNA$QC_Parameter<-rep("nFeature_RNA", times=197)
ggplot(Density__and_box_plot_df_nFeature_RNA, aes(x=QC_Parameter, y=nFeature_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =500)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)

Density__and_box_plot_df_percent_mt<-Density_plot_df[,-c(3,2)]
Density__and_box_plot_df_percent_mt$QC_Parameter<-rep("percent_mt", times=197)
ggplot(Density__and_box_plot_df_percent_mt, aes(x=QC_Parameter, y=percent_mt, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F) + ylim(0, 100)


#Density plots with box plots_log10
#Arrangement of the dataframe
#Density plots
log10_nCount_RNA<-pilot_filt$log10_nCount_RNA
log10_nFeature_RNA<-pilot_filt$log10_nFeature_RNA
Condition<-pilot_filt$Condition
percent_mt<-pilot_filt$percent_mt
Density_plot_df<-data.frame(Condition,log10_nCount_RNA,log10_nFeature_RNA,percent_mt)

Density_plot_df$Cells<-rownames(Density_plot_df)
Density__and_box_plot_df_log10_nCount_RNA<-Density_plot_df[,-c(3,4)]
Density__and_box_plot_df_log10_nCount_RNA$QC_Parameter<-rep("log10_nCount_RNA", times=197)
ggplot(Density__and_box_plot_df_log10_nCount_RNA, aes(x=QC_Parameter, y=log10_nCount_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =2.92942) + geom_hline(yintercept =4.54407)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)


Density__and_box_plot_df_log10_nFeature_RNA<-Density_plot_df[,-c(2,4)]
Density__and_box_plot_df_log10_nFeature_RNA$QC_Parameter<-rep("log10_nFeature_RNA", times=197)
ggplot(Density__and_box_plot_df_log10_nFeature_RNA, aes(x=QC_Parameter, y=log10_nFeature_RNA, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic() + geom_hline(yintercept =2.69897)+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)

Density__and_box_plot_df_percent_mt<-Density_plot_df[,-c(3,2)]
Density__and_box_plot_df_percent_mt$QC_Parameter<-rep("percent_mt", times=197)
ggplot(Density__and_box_plot_df_percent_mt, aes(x=QC_Parameter, y=percent_mt, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .4) +
  geom_boxplot(width = .3,alpha = .4, outlier.shape = NA) +
  scale_color_manual(values=c("#F8766D", "#00BFC4")) + theme_classic()+
  geom_point(aes(color = Condition),
             position = position_jitterdodge(jitter.width = .15,dodge.width = .3),size = 0.8, show.legend = F)+ylim(0,100)














































































#Different quality check graphs

#VIOLIN PLOTS WITH MEAN (nFeature + nCount + percent_mt)
#No filter
VlnPlot(pilot, features= "nFeature_RNA", sort = "decreasing")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5) + scale_y_continuous(limits = c(0,7000))
VlnPlot(pilot, features = c("nCount_RNA"), sort = "decreasing")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5) + scale_y_continuous(limits = c(0,46000))
VlnPlot(pilot, "percent_mt", cols=c("#00BFC4","#F8766D"))+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5) + scale_y_continuous(limits = c(0,100))
VlnPlot(pilot, "percent_cyto",cols=c("#00BFC4","#F8766D"))+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5) + scale_y_continuous(limits = c(0,100))


#Filtered (fixed axis)
VlnPlot(pilot_filt, features= "nFeature_RNA", sort = "Condition")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)+ scale_y_continuous(limits = c(0,7000))+ theme(legend.position = "none")
  #+stat_summary(fun.y = mean, geom='point', colour = "red", width=0.5)
VlnPlot(pilot_filt, features = c("nCount_RNA"), sort = "Condition")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)+ scale_y_continuous(limits = c(0,46000))+ theme(legend.position = "none")
  #+stat_summary(fun.y = mean, geom='point', colour = "red", width=0.5)
-VlnPlot(pilot_filt, features = c("percent_mt"), sort = "Condition")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)+ scale_y_continuous(limits = c(0,100))+ theme(legend.position = "none")
  #+stat_summary(fun.y = mean, geom='point', colour = "red", width=0.5)
VlnPlot(pilot_filt, features = c("percent_cyto"), sort = "decreasing")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5) + scale_y_continuous(limits = c(0,100))+ theme(legend.position = "none")




#Filtered (adapted axis)
P1<-VlnPlot(pilot_filt, features= "nFeature_RNA", sort = "Condition")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)+theme(legend.position = "none")

P2<-VlnPlot(pilot_filt, features = c("nCount_RNA"), sort = "Condition")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)+ theme(legend.position = "none")
P3<-VlnPlot(pilot_filt, features = c("percent_mt"), sort = "Condition")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)+ theme(legend.position = "none")
P4<-VlnPlot(pilot_filt, features = c("percent_cyto"), sort = "decreasing")+
  stat_summary(fun.y = median, geom='crossbar', colour = "black", width=0.5)+ theme(legend.position = "none")

ggarrange(P2,P1,P4,ncol = 3, nrow = 1)

#Statistics


              #Means before filtering
#Fresh cells
Fresh_cells_nFeature_mean<-mean(Fresh_cells$nFeature_RNA, na.rm=TRUE)
Fresh_cells_nFeature_mean

Fresh_cells_nCount_mean<-mean(Fresh_cells$nCount_RNA, na.rm=TRUE)
Fresh_cells_nCount_mean

Fresh_cells_percent_mt_mean<-mean(Fresh_cells$percent_mt, na.rm=TRUE)
Fresh_cells_percent_mt_mean

Fresh_cells_percent_cyto_mean<-mean(Fresh_cells$percent_cyto, na.rm=TRUE)
Fresh_cells_percent_cyto_mean
#Fixed cells
Fixed_cells_nFeature_mean<-mean(Fixed_cells$nFeature_RNA, na.rm=TRUE)
Fixed_cells_nFeature_mean

Fixed_cells_nCount_mean<-mean(Fixed_cells$nCount_RNA, na.rm=TRUE)
Fixed_cells_nCount_mean

Fixed_cells_percent_mt_mean<-mean(Fixed_cells$percent_mt, na.rm=TRUE)
Fixed_cells_percent_mt_mean

Fixed_cells_percent_cyto_mean<-mean(Fixed_cells$percent_cyto, na.rm=TRUE)
Fixed_cells_percent_cyto_mean

              #Means after filtering
#Fresh cells
Fresh_cells_filt_nFeature_mean<-mean(Fresh_cells_filt$nFeature_RNA, na.rm=TRUE)
Fresh_cells_filt_nFeature_mean

Fresh_cells_filt_nCount_mean<-mean(Fresh_cells_filt$nCount_RNA, na.rm=TRUE)
Fresh_cells_filt_nCount_mean

Fresh_cells_filt_percent_mt_mean<-mean(Fresh_cells_filt$percent_mt, na.rm=TRUE)
Fresh_cells_filt_percent_mt_mean

Fresh_cells_filt_percent_cyto_mean<-mean(Fresh_cells_filt$percent_cyto, na.rm=TRUE)
Fresh_cells_filt_percent_cyto_mean
#Fixed cells
Fixed_cells_filt_nFeature_mean<-mean(Fixed_cells_filt$nFeature_RNA, na.rm=TRUE)
Fixed_cells_filt_nFeature_mean

Fixed_cells_filt_nCount_mean<-mean(Fixed_cells_filt$nCount_RNA, na.rm=TRUE)
Fixed_cells_filt_nCount_mean

Fixed_cells_filt_percent_mt_mean<-mean(Fixed_cells_filt$percent_mt, na.rm=TRUE)
Fixed_cells_filt_percent_mt_mean

Fixed_cells_filt_percent_cyto_mean<-mean(Fixed_cells_filt$percent_cyto, na.rm=TRUE)
Fixed_cells_filt_percent_cyto_mean

              #standard deviation before filtering
#Fresh cells
Fresh_cells_nFeature_sd<-sd(Fresh_cells$nFeature_RNA, na.rm=TRUE)
Fresh_cells_nFeature_sd

Fresh_cells_nCount_sd<-sd(Fresh_cells$nCount_RNA, na.rm=TRUE)
Fresh_cells_nCount_sd

Fresh_cells_percent_mt_sd<-sd(Fresh_cells$percent_mt, na.rm=TRUE)
Fresh_cells_percent_mt_sd

Fresh_cells_percent_cyto_sd<-sd(Fresh_cells$percent_cyto, na.rm=TRUE)
Fresh_cells_percent_cyto_sd
#Fixed cells
Fixed_cells_nFeature_sd<-sd(Fixed_cells$nFeature_RNA, na.rm=TRUE)
Fixed_cells_nFeature_sd

Fixed_cells_nCount_sd<-sd(Fixed_cells$nCount_RNA, na.rm=TRUE)
Fixed_cells_nCount_sd

Fixed_cells_percent_mt_sd<-sd(Fixed_cells$percent_mt, na.rm=TRUE)
Fixed_cells_percent_mt_sd

Fixed_cells_percent_cyto_sd<-sd(Fixed_cells$percent_cyto, na.rm=TRUE)
Fixed_cells_percent_cyto_sd

              #standard deviation after filtering
#Fresh cells
Fresh_cells_filt_nFeature_sd<-sd(Fresh_cells_filt$nFeature_RNA, na.rm=TRUE)
Fresh_cells_filt_nFeature_sd

Fresh_cells_filt_nCount_sd<-sd(Fresh_cells_filt$nCount_RNA, na.rm=TRUE)
Fresh_cells_filt_nCount_sd

Fresh_cells_filt_percent_mt_sd<-sd(Fresh_cells_filt$percent_mt, na.rm=TRUE)
Fresh_cells_filt_percent_mt_sd

Fresh_cells_filt_percent_cyto_sd<-sd(Fresh_cells_filt$percent_cyto, na.rm=TRUE)
Fresh_cells_filt_percent_cyto_sd
#Fixed cells
Fixed_cells_filt_nFeature_sd<-sd(Fixed_cells_filt$nFeature_RNA, na.rm=TRUE)
Fixed_cells_filt_nFeature_sd

Fixed_cells_filt_nCount_sd<-sd(Fixed_cells_filt$nCount_RNA, na.rm=TRUE)
Fixed_cells_filt_nCount_sd

Fixed_cells_filt_percent_mt_sd<-sd(Fixed_cells_filt$percent_mt, na.rm=TRUE)
Fixed_cells_filt_percent_mt_sd

Fixed_cells_filt_percent_cyto_sd<-sd(Fixed_cells_filt$percent_cyto, na.rm=TRUE)
Fixed_cells_filt_percent_cyto_sd

#standard error of the mean (SEM)

std <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

            #Fresh cells before filtering
Fresh_cells_nFeature_std<-std(Fresh_cells$nFeature_RNA)
Fresh_cells_nFeature_std

Fresh_cells_nCount_std<-std(Fresh_cells$nCount_RNA, na.rm=TRUE)
Fresh_cells_nCount_std

Fresh_cells_percent_mt_std<-std(Fresh_cells$percent_mt, na.rm=TRUE)
Fresh_cells_percent_mt_std

Fresh_cells_percent_cyto_std<-std(Fresh_cells$percent_cyto, na.rm=TRUE)
Fresh_cells_percent_cyto_std
            #Fixed cells before filtering
Fixed_cells_nFeature_std<-std(Fixed_cells$nFeature_RNA, na.rm=TRUE)
Fixed_cells_nFeature_std

Fixed_cells_nCount_std<-std(Fixed_cells$nCount_RNA, na.rm=TRUE)
Fixed_cells_nCount_std

Fixed_cells_percent_mt_std<-std(Fixed_cells$percent_mt, na.rm=TRUE)
Fixed_cells_percent_mt_std

Fixed_cells_percent_cyto_std<-std(Fixed_cells$percent_cyto, na.rm=TRUE)
Fixed_cells_percent_cyto_std

            #Fresh cells after filtering
Fresh_cells_filt_nFeature_std<-std(Fresh_cells_filt$nFeature_RNA, na.rm=TRUE)
Fresh_cells_filt_nFeature_std

Fresh_cells_filt_nCount_std<-std(Fresh_cells_filt$nCount_RNA, na.rm=TRUE)
Fresh_cells_filt_nCount_std

Fresh_cells_filt_percent_mt_std<-std(Fresh_cells_filt$percent_mt, na.rm=TRUE)
Fresh_cells_filt_percent_mt_std

Fresh_cells_filt_percent_cyto_std<-std(Fresh_cells_filt$percent_cyto, na.rm=TRUE)
Fresh_cells_filt_percent_cyto_std
            #Fixed cells after filtering
Fixed_cells_filt_nFeature_std<-std(Fixed_cells_filt$nFeature_RNA, na.rm=TRUE)
Fixed_cells_filt_nFeature_std

Fixed_cells_filt_nCount_std<-std(Fixed_cells_filt$nCount_RNA, na.rm=TRUE)
Fixed_cells_filt_nCount_std

Fixed_cells_filt_percent_mt_std<-std(Fixed_cells_filt$percent_mt, na.rm=TRUE)
Fixed_cells_filt_percent_mt_std

Fixed_cells_filt_percent_cyto_std<-std(Fixed_cells_filt$percent_cyto, na.rm=TRUE)
Fixed_cells_filt_percent_cyto_std





              #Medians before filtering
#Fresh cells
Fresh_cells_nFeature_median<-median(Fresh_cells$nFeature_RNA, na.rm=TRUE)
Fresh_cells_nFeature_median

Fresh_cells_nCount_median<-median(Fresh_cells$nCount_RNA, na.rm=TRUE)
Fresh_cells_nCount_median

Fresh_cells_percent_mt_median<-median(Fresh_cells$percent_mt, na.rm=TRUE)
Fresh_cells_percent_mt_median

Fresh_cells_percent_cyto_median<-median(Fresh_cells$percent_cyto, na.rm=TRUE)
Fresh_cells_percent_cyto_median
#Fixed cells
Fixed_cells_nFeature_median<-median(Fixed_cells$nFeature_RNA, na.rm=TRUE)
Fixed_cells_nFeature_median

Fixed_cells_nCount_median<-median(Fixed_cells$nCount_RNA, na.rm=TRUE)
Fixed_cells_nCount_median

Fixed_cells_percent_mt_median<-median(Fixed_cells$percent_mt, na.rm=TRUE)
Fixed_cells_percent_mt_median

Fixed_cells_percent_cyto_median<-median(Fixed_cells$percent_cyto, na.rm=TRUE)
Fixed_cells_percent_cyto_median

              #Medians after filtering
#Fresh cells
Fresh_cells_filt_nFeature_median<-median(Fresh_cells_filt$nFeature_RNA, na.rm=TRUE)
Fresh_cells_filt_nFeature_median

Fresh_cells_filt_nCount_median<-median(Fresh_cells_filt$nCount_RNA, na.rm=TRUE)
Fresh_cells_filt_nCount_median

Fresh_cells_filt_percent_mt_median<-median(Fresh_cells_filt$percent_mt, na.rm=TRUE)
Fresh_cells_filt_percent_mt_median

Fresh_cells_filt_percent_cyto_median<-median(Fresh_cells_filt$percent_cyto, na.rm=TRUE)
Fresh_cells_filt_percent_cyto_median
#Fixed cells
Fixed_cells_filt_nFeature_median<-median(Fixed_cells_filt$nFeature_RNA, na.rm=TRUE)
Fixed_cells_filt_nFeature_median

Fixed_cells_filt_nCount_median<-median(Fixed_cells_filt$nCount_RNA, na.rm=TRUE)
Fixed_cells_filt_nCount_median

Fixed_cells_filt_percent_mt_median<-median(Fixed_cells_filt$percent_mt, na.rm=TRUE)
Fixed_cells_filt_percent_mt_median

Fixed_cells_filt_percent_cyto_median<-median(Fixed_cells_filt$percent_cyto, na.rm=TRUE)
Fixed_cells_filt_percent_cyto_median









#Exporting the data to test it on Sigma Plot


      #Fresh cells-           Before filtering

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fresh_cells$nCount_RNA, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fresh_cells$nCount_RNA", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fresh_cells$nFeature_RNA, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fresh_cells$nFeature_RNA", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fresh_cells$percent_mt, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fresh_cells$percent_mt", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fresh_cells$percent_cyto, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fresh_cells$percent_cyto", col.names = TRUE, row.names = TRUE, append = FALSE)

      #Fixed cells-           Before filtering

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fixed_cells$nCount_RNA, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fixed_cells$nCount_RNA", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fixed_cells$nFeature_RNA, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fixed_cells$nFeature_RNA", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fixed_cells$percent_mt, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fixed_cells$percent_mt", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fixed_cells$percent_cyto, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fixed_cells$percent_cyto", col.names = TRUE, row.names = TRUE, append = FALSE)

      #Fresh cells-           After filtering

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fresh_cells_filt$nCount_RNA, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fresh_cells_filt$nCount_RNA", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fresh_cells_filt$nFeature_RNA, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fresh_cells_filt$nFeature_RNA", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fresh_cells_filt$percent_mt, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fresh_cells_filt$percent_mt", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fresh_cells_filt$percent_cyto, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fresh_cells_filt$percent_cyto", col.names = TRUE, row.names = TRUE, append = FALSE)

      #Fixed cells_filt-           After filtering

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fixed_cells_filt$nCount_RNA, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fixed_cells_filt$nCount_RNA", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fixed_cells_filt$nFeature_RNA, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fixed_cells_filt$nFeature_RNA", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fixed_cells_filt$percent_mt, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fixed_cells_filt$percent_mt", col.names = TRUE, row.names = TRUE, append = FALSE)

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Fixed_cells_filt$percent_cyto, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/Parameters exported from R/", paste(time.name), "Parameters.xlsx"))
           , sheetName = "Fixed_cells_filt$percent_cyto", col.names = TRUE, row.names = TRUE, append = FALSE)





#Normality test
#The null hypothesis of these tests is that “sample distribution is normal”. If the test is significant, the distribution is non-normal.
#Before filtering (BF)
    #Fresh cells
ggdensity(Fresh_cells$nFeature_RNA, 
          main = "Density plot of Before filetring Fresh cells nFeature",
          xlab = "nFeature")
shapiro.test(Fresh_cells$nFeature_RNA)

ggdensity(Fresh_cells$nCount_RNA, 
          main = "Density plot of Before filetring Fresh cells nCount_RNA",
          xlab = "nCount_RNA")
shapiro.test(Fresh_cells$nCount_RNA)

ggdensity(Fresh_cells$percent_mt, 
          main = "Density plot of Before filetring Fresh cells percent_mt",
          xlab = "percent_mt")
shapiro.test(Fresh_cells$percent_mt)

ggdensity(Fresh_cells$percent_cyto, 
          main = "Density plot of Before filetring Fresh cells percent_cyto",
          xlab = "percent_cyto")
shapiro.test(Fresh_cells$percent_cyto)

    #Fixed cells

ggdensity(Fixed_cells$nFeature_RNA, 
          main = "Density plot of Before filetring Fixed cells nFeature",
          xlab = "nFeature")
shapiro.test(Fixed_cells$nFeature_RNA)

ggdensity(Fixed_cells$nCount_RNA, 
          main = "Density plot of Before filetring Fixed cells nCount_RNA",
          xlab = "nCount_RNA")
shapiro.test(Fixed_cells$nCount_RNA)

ggdensity(Fixed_cells$percent_mt, 
          main = "Density plot of Before filetring Fixed cells percent_mt",
          xlab = "percent_mt")
shapiro.test(Fixed_cells$percent_mt)

ggdensity(Fixed_cells$percent_cyto, 
          main = "Density plot of Before filetring Fixed cells percent_cyto",
          xlab = "percent_cyto")
shapiro.test(Fixed_cells$percent_cyto)

#After filtering 

    #Fresh cells
ggdensity(Fresh_cells_filt$nFeature_RNA, 
          main = "Density plot of After filetring Fresh cells nFeature",
          xlab = "nFeature")
shapiro.test(Fresh_cells_filt$nFeature_RNA)

ggdensity(Fresh_cells_filt$nCount_RNA, 
          main = "Density plot of After filetring Fresh cells nCount_RNA",
          xlab = "nCount_RNA")
shapiro.test(Fresh_cells_filt$nCount_RNA)

ggdensity(Fresh_cells_filt$percent_mt, 
          main = "Density plot of After filetring Fresh cells percent_mt",
          xlab = "percent_mt")
shapiro.test(Fresh_cells_filt$percent_mt)

ggdensity(Fresh_cells_filt$percent_cyto, 
          main = "Density plot of After filetring Fresh cells percent_cyto",
          xlab = "percent_cyto")
shapiro.test(Fresh_cells_filt$percent_cyto)

    #Fixed cells

ggdensity(Fixed_cells_filt$nFeature_RNA, 
          main = "Density plot of After filetring Fixed cells nFeature",
          xlab = "nFeature")
shapiro.test(Fixed_cells_filt$nFeature_RNA)

ggdensity(Fixed_cells_filt$nCount_RNA, 
          main = "Density plot of After filetring Fixed cells nCount_RNA",
          xlab = "nCount_RNA")
shapiro.test(Fixed_cells_filt$nCount_RNA)

ggdensity(Fixed_cells_filt$percent_mt, 
          main = "Density plot of After filetring Fixed cells percent_mt",
          xlab = "percent_mt")
shapiro.test(Fixed_cells_filt$percent_mt)

ggdensity(Fixed_cells_filt$percent_cyto, 
          main = "Density plot of After filetring Fixed cells percent_cyto",
          xlab = "percent_cyto")
shapiro.test(Fixed_cells_filt$percent_cyto)

#Wilcoxon rank-Sum test

      #Fresh Cells Before - Fixed cells Before
wilcox.test(Fresh_cells$nFeature_RNA, Fixed_cells$nFeature_RNA, paired=FALSE)
wilcox.test(Fresh_cells$nCount_RNA, Fixed_cells$nCount_RNA, paired=FALSE)
wilcox.test(Fresh_cells$percent_mt, Fixed_cells$percent_mt, paired=FALSE)
wilcox.test(Fresh_cells$percent_cyto, Fixed_cells$percent_cyto, paired=FALSE)

      #Fresh Cells After - Fixed cells After
wilcox.test(Fresh_cells_filt$nFeature_RNA, Fixed_cells_filt$nFeature_RNA, paired=FALSE)
wilcox.test(Fresh_cells_filt$nCount_RNA, Fixed_cells_filt$nCount_RNA, paired=FALSE)
wilcox.test(Fresh_cells_filt$percent_mt, Fixed_cells_filt$percent_mt, paired=FALSE)
wilcox.test(Fresh_cells_filt$percent_cyto, Fixed_cells_filt$percent_cyto, paired=FALSE)


#Dot plot mitochondrial genes expression
mito_genes_to_plot <- grep("^mt",pilot@assays[["RNA"]]@data@Dimnames[[1]], value=TRUE) 
DotPlot(object= pilot, features = mito_genes_to_plot, dot.scale = 6)+
  RotatedAxis()
DotPlot(object= pilot_filt, features = mito_genes_to_plot, dot.scale = 5,scale=FALSE)+
  RotatedAxis()











#To check the real differential expressed mitochondrial genes between fresh and fixed
#I do a reduced dataset with only yhe mitochondrial genes to do the DE analysis 
#and obtain which ones are differential expressed in fresh
pilot.filt.df<-data.frame(pilot_filt@assays[["RNA"]]@counts)
mito2 <- grep("^mt",rownames(pilot.filt.df))
Pilot.mito<-pilot.filt.df[mito2,]
Pilot_mito<-CreateSeuratObject(counts = Pilot.mito, project = "Pilot analysis", min.cells = 0, min.features = 0)

current.ids2 <- c("JEN.AU.s004","JEN.AU.s002")
new.ids2 <- c("Fresh_cells","Fixed_cells")
Pilot_mito@active.ident <- plyr::mapvalues(x = Pilot_mito@active.ident, from = current.ids2, to = new.ids2)
Pilot_mito$Condition <- Idents(Pilot_mito)

Pilot_mito$Condition <- factor(x=Pilot_mito$Condition, levels=c("Fresh_cells","Fixed_cells"))




Mito_DE<-FindMarkers(Pilot_mito,ident.1 ="Fresh_cells",ident.2 ="Fixed_cells",assay="RNA", slot="counts")
Mito_DE
Mito_DE$BH_p_val_adj<-p.adjust(Mito_DE$p_val,method="BH")
Mito_DE



#Now I plot with the assay=RNA
DotPlot(object= Pilot_mito, features = mito_genes_to_plot,scale=FALSE, dot.scale = 5, assay= "RNA")+
  RotatedAxis()


Mito_DE %>%
  slice_max(n =10, order_by = avg_log2FC) -> top10_log2FC
top10_log2FC_HM<-DoHeatmap(Pilot_mito, features = rownames(top10_log2FC), label=FALSE,slot="counts")
top10_log2FC_HM

DoHeatmap(object= Pilot_mito, features = mito_genes_to_plot, slot="counts")

#here checking some means to see if the results of the DE analysis 
#or the ones from the graph are the correct ones
Fresh.df<-data.frame(Fresh_cells_filt@assays[["RNA"]]@counts)
Fixed.df<-data.frame(Fixed_cells_filt@assays[["RNA"]]@counts)


Rnr2<- grep("mt-Rnr2",rownames(Fresh.df))
Fresh_mt_Rnr2<-Matrix::colSums(Fresh.df[Rnr2,])
mean(Fresh_mt_Rnr2)
Rnr2_<- grep("mt-Rnr2",rownames(Fixed.df))
Fixed_mt_Rnr2<-Matrix::colSums(Fixed.df[Rnr2_,])
mean(Fixed_mt_Rnr2)


mtNd5<- grep("mtNd5",rownames(Fresh.df))
Fresh_mt_mtNd5<-Matrix::colSums(Fresh.df[mtNd5,])
mean(Fresh_mt_mtNd5)
mtNd5_<- grep("mtNd5",rownames(Fixed.df))
Fixed_mt_mtNd5<-Matrix::colSums(Fixed.df[mtNd5_,])
mean(Fixed_mt_mtNd5)

mtRnr1<- grep("mt-Rnr1",rownames(Fresh.df))
Fresh_mt_mtRnr1<-Matrix::colSums(Fresh.df[mtRnr1,])
mean(Fresh_mt_mtRnr1)
mtRnr1_<- grep("mt-Rnr1",rownames(Fixed.df))
Fixed_mt_mtRnr1<-Matrix::colSums(Fixed.df[mtRnr1_,])
mean(Fixed_mt_mtRnr1)

#Dot plot to show gene expression

mito_genes <- c("mt-Rnr1", "mt-Rnr2", "mt-Ta",   "mt-Tc" ,  "mt-Td" ,  "mt-Te" ,  "mt-Tf" ,  "mt-Tg"  , "mt-Th" ,  "mt-Ti" ,  "mt-Tk" ,
                  "mt-Tl1" , "mt-Tl2" , "mt-Tm" ,  "mt-Tn" , "mt-Tp" ,  "mt-Tq" ,  "mt-Tr" ,  "mt-Ts1",  "mt-Ts2" , "mt-Tt" ,  "mt-Tv" ,
                  "mt-Tw" ,  "mt-Ty" ,  "mtAtp6" , "mtAtp8" , "mtCo1" ,  "mtCo2"  , "mtCo3" ,  "mtCytb" ,
                  "mtNd1" ,  "mtNd2" ,  "mtNd3" ,  "mtNd4" ,  "mtNd4l" , "mtNd5" ,  "mtNd6"  )
dittoDotPlot(Pilot_mito, mito_genes,group.by="Condition",assay="RNA", slot="counts", scale=FALSE,y.reorder=c(2,1),
             max.color =  "#D73027", min.color= "#D0D0D1")


#Export dataframe to test on Sigma Plot

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Pilot.mito, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Pilot.mito.xlsx"))
           , sheetName = "Pilot.mito", col.names = TRUE, row.names = TRUE, append = FALSE)



#To obtain the full list of genes from FindMarkers
options(scipen = 999)
Mito_DE_all<-FindMarkers(Pilot_mito,ident.1 ="Fresh_cells",ident.2 ="Fixed_cells",assay="RNA", slot="counts",  
                     min.cells.group = 1, 
                     min.cells.feature = 1,
                     min.pct = 0,
                     logfc.threshold = 0,
                     only.pos = FALSE)
Mito_DE_all
Mito_DE_all$BH_p_val_adj<-p.adjust(Mito_DE_all$p_val,method="BH")
Mito_DE_all


time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Mito_DE_all, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Mito_DE.xlsx"))
           , sheetName = "Mito_DE", col.names = TRUE, row.names = TRUE, append = FALSE)





































#Checking stress genes (Hspb1,Egr1 and Bax)


VlnPlot(pilot, features = c("Hspb1", "Egr1", "Bax"),sort = "decreasing")
VlnPlot(pilot_filt, features = c("Hspb1", "Egr1", "Bax"), sort = "Condition")

DotPlot(pilot, features = c("Hspb1", "Egr1", "Bax"))
DotPlot(pilot_filt, features = c("Hspb1", "Egr1", "Bax"))

      #Statistics
#first I create an columns on metadata with only the expression of Hspb1, Egr1 and Bax

#Pilot_stress<-pilot[(pilot@assays[["RNA"]]@counts@Dimnames[[1]]=="Hspb1" | pilot@assays[["RNA"]]@counts@Dimnames[[1]]=="Egr1"|pilot@assays[["RNA"]]@counts@Dimnames[[1]]=="Bax"),]
#Pilot_filt_stress<-pilot_filt[(pilot_filt@assays[["RNA"]]@counts@Dimnames[[1]]=="Hspb1" | pilot_filt@assays[["RNA"]]@counts@Dimnames[[1]]=="Egr1"|pilot_filt@assays[["RNA"]]@counts@Dimnames[[1]]=="Bax"),]

pilot[["Hspb1"]] <- FetchData(pilot, vars="Hspb1")
pilot[["Egr1"]] <- FetchData(pilot, vars="Egr1")
pilot[["Bax"]] <- FetchData(pilot, vars="Bax")

pilot_filt[["Hspb1"]] <- FetchData(pilot_filt, vars="Hspb1")
pilot_filt[["Egr1"]] <- FetchData(pilot_filt, vars="Egr1")
pilot_filt[["Bax"]] <- FetchData(pilot_filt, vars="Bax")

#Normality tests

#Before filtering
ggdensity(pilot$Hspb1, 
          main = "Density plot of Before filetring pilot Hspb1",
          xlab = "Level of expression")
shapiro.test(pilot$Hspb1)

ggdensity(pilot$Egr1, 
          main = "Density plot of Before filetring pilot Egr1",
          xlab = "Level of expression")
shapiro.test(pilot$Egr1)

ggdensity(pilot$Bax, 
          main = "Density plot of Before filetring pilot Bax",
          xlab = "Level of expression")
shapiro.test(pilot$Bax)


#After filtering
ggdensity(pilot_filt$Hspb1, 
          main = "Density plot of Before filetring pilot_filt Hspb1",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Hspb1)

ggdensity(pilot_filt$Egr1, 
          main = "Density plot of Before filetring pilot_filt Egr1",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Egr1)

ggdensity(pilot_filt$Bax, 
          main = "Density plot of Before filetring pilot_filt Bax",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Bax)


#Wilcoxon rank-Sum test
#To compare I need to split
#Before
pilot.list<-SplitObject(pilot, split.by = "Condition")

Fresh_cells<-pilot.list$Fresh_cells
Fixed_cells<-pilot.list$Fixed_cells
#After
pilot.list_filt<-SplitObject(pilot_filt, split.by = "Condition")

Fresh_cells_filt<-pilot.list_filt$Fresh_cells
Fixed_cells_filt<-pilot.list_filt$Fixed_cells

#Fresh Cells Before - Fixed cells Before
wilcox.test(Fresh_cells$Hspb1, Fixed_cells$Hspb1, paired=FALSE)
wilcox.test(Fresh_cells$Egr1, Fixed_cells$Egr1, paired=FALSE)
wilcox.test(Fresh_cells$Bax, Fixed_cells$Bax, paired=FALSE)

#Fresh Cells After - Fixed cells After
wilcox.test(Fresh_cells_filt$Hspb1, Fixed_cells_filt$Hspb1, paired=FALSE)
wilcox.test(Fresh_cells_filt$Egr1, Fixed_cells_filt$Egr1, paired=FALSE)
wilcox.test(Fresh_cells_filt$Bax, Fixed_cells_filt$Bax, paired=FALSE)


#Means before filtering
#Fresh cells
Fresh_cells_Hspb1_mean<-mean(Fresh_cells$Hspb1, na.rm=TRUE)
Fresh_cells_Hspb1_mean
Fresh_cells_Egr1_mean<-mean(Fresh_cells$Egr1, na.rm=TRUE)
Fresh_cells_Egr1_mean
Fresh_cells_Bax_mean<-mean(Fresh_cells$Bax, na.rm=TRUE)
Fresh_cells_Bax_mean

#Fixed cells
Fixed_cells_Hspb1_mean<-mean(Fixed_cells$Hspb1, na.rm=TRUE)
Fixed_cells_Hspb1_mean
Fixed_cells_Egr1_mean<-mean(Fixed_cells$Egr1, na.rm=TRUE)
Fixed_cells_Egr1_mean
Fixed_cells_Bax_mean<-mean(Fixed_cells$Bax, na.rm=TRUE)
Fixed_cells_Bax_mean

#Means after filtering
#Fresh cells
Fresh_cells_filt_Hspb1_mean<-mean(Fresh_cells_filt$Hspb1, na.rm=TRUE)
Fresh_cells_filt_Hspb1_mean
Fresh_cells_filt_Egr1_mean<-mean(Fresh_cells_filt$Egr1, na.rm=TRUE)
Fresh_cells_filt_Egr1_mean
Fresh_cells_filt_Bax_mean<-mean(Fresh_cells_filt$Bax, na.rm=TRUE)
Fresh_cells_filt_Bax_mean

#Fixed cells
Fixed_cells_filt_Hspb1_mean<-mean(Fixed_cells_filt$Hspb1, na.rm=TRUE)
Fixed_cells_filt_Hspb1_mean
Fixed_cells_filt_Egr1_mean<-mean(Fixed_cells_filt$Egr1, na.rm=TRUE)
Fixed_cells_filt_Egr1_mean
Fixed_cells_filt_Bax_mean<-mean(Fixed_cells_filt$Bax, na.rm=TRUE)
Fixed_cells_filt_Bax_mean

#Summary means
Stress_genes<-c("Hspb1", "Egr1", "Bax")
Fresh_cells_Before_filt_mean<-c(Fresh_cells_Hspb1_mean,Fresh_cells_Egr1_mean,Fresh_cells_Bax_mean)
Fixed_cells_Before_filt_mean<-c(Fixed_cells_Hspb1_mean,Fixed_cells_Egr1_mean,Fixed_cells_Bax_mean)
Fresh_cells_After_filt_mean<-c(Fresh_cells_filt_Hspb1_mean,Fresh_cells_filt_Egr1_mean,Fresh_cells_filt_Bax_mean)
Fixed_cells_After_filt_mean<-c(Fixed_cells_filt_Hspb1_mean,Fixed_cells_filt_Egr1_mean,Fixed_cells_filt_Bax_mean)
Stress_genes_summary_mean<-data.frame(Stress_genes,Fresh_cells_Before_filt_mean,Fixed_cells_Before_filt_mean,Fresh_cells_After_filt_mean,Fixed_cells_After_filt_mean)
Stress_genes_summary_mean

#standard error of the mean (SEM)

std <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

#stds before filtering
#Fresh cells
Fresh_cells_Hspb1_std<-std(Fresh_cells$Hspb1, na.rm=TRUE)
Fresh_cells_Hspb1_std
Fresh_cells_Egr1_std<-std(Fresh_cells$Egr1, na.rm=TRUE)
Fresh_cells_Egr1_std
Fresh_cells_Bax_std<-std(Fresh_cells$Bax, na.rm=TRUE)
Fresh_cells_Bax_std

#Fresh cells
Fixed_cells_Hspb1_std<-std(Fixed_cells$Hspb1, na.rm=TRUE)
Fixed_cells_Hspb1_std
Fixed_cells_Egr1_std<-std(Fixed_cells$Egr1, na.rm=TRUE)
Fixed_cells_Egr1_std
Fixed_cells_Bax_std<-std(Fixed_cells$Bax, na.rm=TRUE)
Fixed_cells_Bax_std

#stds after filtering
#Fresh cells
Fresh_cells_filt_Hspb1_std<-std(Fresh_cells_filt$Hspb1, na.rm=TRUE)
Fresh_cells_filt_Hspb1_std
Fresh_cells_filt_Egr1_std<-std(Fresh_cells_filt$Egr1, na.rm=TRUE)
Fresh_cells_filt_Egr1_std
Fresh_cells_filt_Bax_std<-std(Fresh_cells_filt$Bax, na.rm=TRUE)
Fresh_cells_filt_Bax_std

#Fresh cells
Fixed_cells_filt_Hspb1_std<-std(Fixed_cells_filt$Hspb1, na.rm=TRUE)
Fixed_cells_filt_Hspb1_std
Fixed_cells_filt_Egr1_std<-std(Fixed_cells_filt$Egr1, na.rm=TRUE)
Fixed_cells_filt_Egr1_std
Fixed_cells_filt_Bax_std<-std(Fixed_cells_filt$Bax, na.rm=TRUE)
Fixed_cells_filt_Bax_std

#Summary stds
Stress_genes<-c("Hspb1", "Egr1", "Bax")
Fresh_cells_Before_filt_std<-c(Fresh_cells_Hspb1_std,Fresh_cells_Egr1_std,Fresh_cells_Bax_std)
Fixed_cells_Before_filt_std<-c(Fixed_cells_Hspb1_std,Fixed_cells_Egr1_std,Fixed_cells_Bax_std)
Fresh_cells_After_filt_std<-c(Fresh_cells_filt_Hspb1_std,Fresh_cells_filt_Egr1_std,Fresh_cells_filt_Bax_std)
Fixed_cells_After_filt_std<-c(Fixed_cells_filt_Hspb1_std,Fixed_cells_filt_Egr1_std,Fixed_cells_filt_Bax_std)
Stress_genes_summary_std<-data.frame(Stress_genes,Fresh_cells_Before_filt_std,Fixed_cells_Before_filt_std,Fresh_cells_After_filt_std,Fixed_cells_After_filt_std)
Stress_genes_summary_std

#Summary stress genes statistics

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Stress_genes_summary_mean, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Stress_genes_summary_mean", col.names = TRUE, row.names = TRUE)
write.xlsx(Stress_genes_summary_std, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Stress_genes_summary_std", col.names = TRUE, row.names = TRUE, append = TRUE)

#Summary stress genes expression for box plots
#Before
write.xlsx(Fresh_cells$Hspb1, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fresh_cells_Before filt_Hspb1", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Fixed_cells$Hspb1, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fixed_cells_Before filt_Hspb1", col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx(Fresh_cells$Egr1, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fresh_cells_Before filt_Egr1", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Fixed_cells$Egr1, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fixed_cells_Before filt_Egr1", col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx(Fresh_cells$Bax, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fresh_cells_Before filt_Bax", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Fixed_cells$Bax, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fixed_cells_Before filt_Bax", col.names = TRUE, row.names = TRUE, append = TRUE)

#After

write.xlsx(Fresh_cells_filt$Hspb1, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fresh_cells_After filt_Hspb1", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Fixed_cells_filt$Hspb1, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fixed_cells_After filt_Hspb1", col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx(Fresh_cells_filt$Egr1, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fresh_cells_After filt_Egr1", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Fixed_cells_filt$Egr1, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fixed_cells_After filt_Egr1", col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx(Fresh_cells_filt$Bax, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fresh_cells_After filt_Bax", col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Fixed_cells_filt$Bax, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary.xlsx"))
           , sheetName = "Fixed_cells_After filt_Bax", col.names = TRUE, row.names = TRUE, append = TRUE)




#Increasing stress genes study with ("Atf3", "Junb", "Ier3", "Socs3", "P2ry13", "Bik", "Casp11", "C3", "Cd63", "Cldn1", "Rxrg", "Arc", "Rem2", "Hmbs", "Sdha", "Gapdh", "p53", "casp9")


VlnPlot(pilot, features = c("Atf3", "Junb", "Ier3", "Socs3", "P2ry13", "Bik", "C3", "Cd63", "Cldn1", "Rxrg", 
                            "Arc", "Rem2", "Hmbs", "Sdha", "Gapdh", "p53", "Casp9","Casp1"), sort = "active.ident", ncol=7)
VlnPlot(pilot_filt, features = c("Atf3", "Junb", "Ier3", "Socs3", "P2ry13", "Bik", "C3", "Cd63", "Cldn1", "Rxrg", 
                                 "Arc", "Rem2", "Hmbs", "Sdha", "Gapdh", "p53", "Casp9","Casp1"), sort = "Condition", ncol=7)

gene_list<-row.names(pilot@assays[["RNA"]])

gene_list<-as.matrix(gene_list)


#Statistics
#first I create an columns on metadata with only the expression of the genes

pilot[["Atf3"]] <- FetchData(pilot, vars="Atf3")
pilot[["Junb"]] <- FetchData(pilot, vars="Junb")
pilot[["Ier3"]] <- FetchData(pilot, vars="Ier3")
pilot[["Socs3"]] <- FetchData(pilot, vars="Socs3")
pilot[["P2ry13"]] <- FetchData(pilot, vars="P2ry13")
pilot[["Cd63"]] <- FetchData(pilot, vars="Cd63")
pilot[["Cldn1"]] <- FetchData(pilot, vars="Cldn1")
pilot[["Rxrg"]] <- FetchData(pilot, vars="Rxrg")
pilot[["Arc"]] <- FetchData(pilot, vars="Arc")
pilot[["Rem2"]] <- FetchData(pilot, vars="Rem2")
pilot[["Hmbs"]] <- FetchData(pilot, vars="Hmbs")
pilot[["Sdha"]] <- FetchData(pilot, vars="Sdha")
pilot[["Gapdh"]] <- FetchData(pilot, vars="Gapdh")
pilot[["Casp9"]] <- FetchData(pilot, vars="Casp9")
pilot[["Casp1"]] <- FetchData(pilot, vars="Casp1")

pilot_filt[["Atf3"]] <- FetchData(pilot_filt, vars="Atf3")
pilot_filt[["Junb"]] <- FetchData(pilot_filt, vars="Junb")
pilot_filt[["Ier3"]] <- FetchData(pilot_filt, vars="Ier3")
pilot_filt[["Socs3"]] <- FetchData(pilot_filt, vars="Socs3")
pilot_filt[["P2ry13"]] <- FetchData(pilot_filt, vars="P2ry13")
pilot_filt[["Cd63"]] <- FetchData(pilot_filt, vars="Cd63")
pilot_filt[["Cldn1"]] <- FetchData(pilot_filt, vars="Cldn1")
pilot_filt[["Rxrg"]] <- FetchData(pilot_filt, vars="Rxrg")
pilot_filt[["Arc"]] <- FetchData(pilot_filt, vars="Arc")
pilot_filt[["Rem2"]] <- FetchData(pilot_filt, vars="Rem2")
pilot_filt[["Hmbs"]] <- FetchData(pilot_filt, vars="Hmbs")
pilot_filt[["Sdha"]] <- FetchData(pilot_filt, vars="Sdha")
pilot_filt[["Gapdh"]] <- FetchData(pilot_filt, vars="Gapdh")
pilot_filt[["Casp9"]] <- FetchData(pilot_filt, vars="Casp9")
pilot_filt[["Casp1"]] <- FetchData(pilot_filt, vars="Casp1")


#Normality tests

#Before filtering
ggdensity(pilot$Atf3, 
          main = "Density plot of Before filetring pilot Atf3",
          xlab = "Level of expression")
shapiro.test(pilot$Atf3)

ggdensity(pilot$Junb, 
          main = "Density plot of Before filetring pilot Junb",
          xlab = "Level of expression")
shapiro.test(pilot$Junb)

ggdensity(pilot$Ier3, 
          main = "Density plot of Before filetring pilot Ier3",
          xlab = "Level of expression")
shapiro.test(pilot$Ier3)

ggdensity(pilot$Socs3, 
          main = "Density plot of Before filetring pilot Socs3",
          xlab = "Level of expression")
shapiro.test(pilot$Socs3)

ggdensity(pilot$P2ry13, 
          main = "Density plot of Before filetring pilot P2ry13",
          xlab = "Level of expression")
shapiro.test(pilot$P2ry13)


ggdensity(pilot$Cd63, 
          main = "Density plot of Before filetring pilot Cd63",
          xlab = "Level of expression")
shapiro.test(pilot$Cd63)

ggdensity(pilot$Cldn1, 
          main = "Density plot of Before filetring pilot Cldn1",
          xlab = "Level of expression")
shapiro.test(pilot$Cldn1)

ggdensity(pilot$Rxrg, 
          main = "Density plot of Before filetring pilot Rxrg",
          xlab = "Level of expression")
shapiro.test(pilot$Rxrg)

ggdensity(pilot$Arc, 
          main = "Density plot of Before filetring pilot Arc",
          xlab = "Level of expression")
shapiro.test(pilot$Arc)

ggdensity(pilot$Rem2, 
          main = "Density plot of Before filetring pilot Rem2",
          xlab = "Level of expression")
shapiro.test(pilot$Rem2)

ggdensity(pilot$Hmbs, 
          main = "Density plot of Before filetring pilot Hmbs",
          xlab = "Level of expression")
shapiro.test(pilot$Hmbs)

ggdensity(pilot$Sdha, 
          main = "Density plot of Before filetring pilot Sdha",
          xlab = "Level of expression")
shapiro.test(pilot$Sdha)

ggdensity(pilot$Gapdh, 
          main = "Density plot of Before filetring pilot Gapdh",
          xlab = "Level of expression")
shapiro.test(pilot$Gapdh)

ggdensity(pilot$Casp9, 
          main = "Density plot of Before filetring pilot Casp9",
          xlab = "Level of expression")
shapiro.test(pilot$Casp9)
ggdensity(pilot$Casp1, 
          main = "Density plot of Before filetring pilot Casp1",
          xlab = "Level of expression")
shapiro.test(pilot$Casp1)


#After filtering
ggdensity(pilot_filt$Atf3, 
          main = "Density plot of Before filetring pilot_filt Atf3",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Atf3)

ggdensity(pilot_filt$Junb, 
          main = "Density plot of Before filetring pilot_filt Junb",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Junb)

ggdensity(pilot_filt$Ier3, 
          main = "Density plot of Before filetring pilot_filt Ier3",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Ier3)

ggdensity(pilot_filt$Socs3, 
          main = "Density plot of Before filetring pilot_filt Socs3",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Socs3)

ggdensity(pilot_filt$P2ry13, 
          main = "Density plot of Before filetring pilot_filt P2ry13",
          xlab = "Level of expression")
shapiro.test(pilot_filt$P2ry13)


ggdensity(pilot_filt$Cd63, 
          main = "Density plot of Before filetring pilot_filt Cd63",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Cd63)

ggdensity(pilot_filt$Cldn1, 
          main = "Density plot of Before filetring pilot_filt Cldn1",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Cldn1)

ggdensity(pilot_filt$Rxrg, 
          main = "Density plot of Before filetring pilot_filt Rxrg",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Rxrg)

ggdensity(pilot_filt$Arc, 
          main = "Density plot of Before filetring pilot_filt Arc",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Arc)

ggdensity(pilot_filt$Rem2, 
          main = "Density plot of Before filetring pilot_filt Rem2",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Rem2)

ggdensity(pilot_filt$Hmbs, 
          main = "Density plot of Before filetring pilot_filt Hmbs",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Hmbs)

ggdensity(pilot_filt$Sdha, 
          main = "Density plot of Before filetring pilot_filt Sdha",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Sdha)

ggdensity(pilot_filt$Gapdh, 
          main = "Density plot of Before filetring pilot_filt Gapdh",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Gapdh)

ggdensity(pilot_filt$Casp9, 
          main = "Density plot of Before filetring pilot_filt Casp9",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Casp9)
ggdensity(pilot_filt$Casp1, 
          main = "Density plot of Before filetring pilot_filt Casp1",
          xlab = "Level of expression")
shapiro.test(pilot_filt$Casp1)


#Wilcoxon rank-Sum test
#To compare I need to split
#Before
pilot.list<-SplitObject(pilot, split.by = "Condition")

Fresh_cells<-pilot.list$Fresh_cells
Fixed_cells<-pilot.list$Fixed_cells
#After
pilot.list_filt<-SplitObject(pilot_filt, split.by = "Condition")

Fresh_cells_filt<-pilot.list_filt$Fresh_cells
Fixed_cells_filt<-pilot.list_filt$Fixed_cells

#Fresh Cells Before - Fixed cells Before
wilcox.test(Fresh_cells$Atf3, Fixed_cells$Atf3, paired=FALSE)
wilcox.test(Fresh_cells$Junb, Fixed_cells$Junb, paired=FALSE)
wilcox.test(Fresh_cells$Ier3, Fixed_cells$Ier3, paired=FALSE)
wilcox.test(Fresh_cells$Socs3, Fixed_cells$Socs3, paired=FALSE)
wilcox.test(Fresh_cells$P2ry13, Fixed_cells$P2ry13, paired=FALSE)
wilcox.test(Fresh_cells$Cd63, Fixed_cells$Cd63, paired=FALSE)
wilcox.test(Fresh_cells$Cldn1, Fixed_cells$Cldn1, paired=FALSE)
wilcox.test(Fresh_cells$Rxrg, Fixed_cells$Rxrg, paired=FALSE)
wilcox.test(Fresh_cells$Arc, Fixed_cells$Arc, paired=FALSE)
wilcox.test(Fresh_cells$Rem2, Fixed_cells$Rem2, paired=FALSE)
wilcox.test(Fresh_cells$Hmbs, Fixed_cells$Hmbs, paired=FALSE)
wilcox.test(Fresh_cells$Sdha, Fixed_cells$Sdha, paired=FALSE)
wilcox.test(Fresh_cells$Gapdh, Fixed_cells$Gapdh, paired=FALSE)
wilcox.test(Fresh_cells$Casp9, Fixed_cells$Casp9, paired=FALSE)
wilcox.test(Fresh_cells$Casp1, Fixed_cells$Casp1, paired=FALSE)

#Fresh Cells After - Fixed cells After
wilcox.test(Fresh_cells_filt$Atf3, Fixed_cells_filt$Atf3, paired=FALSE)
wilcox.test(Fresh_cells_filt$Junb, Fixed_cells_filt$Junb, paired=FALSE)
wilcox.test(Fresh_cells_filt$Ier3, Fixed_cells_filt$Ier3, paired=FALSE)
wilcox.test(Fresh_cells_filt$Socs3, Fixed_cells_filt$Socs3, paired=FALSE)
wilcox.test(Fresh_cells_filt$P2ry13, Fixed_cells_filt$P2ry13, paired=FALSE)
wilcox.test(Fresh_cells_filt$Cd63, Fixed_cells_filt$Cd63, paired=FALSE)
wilcox.test(Fresh_cells_filt$Cldn1, Fixed_cells_filt$Cldn1, paired=FALSE)
wilcox.test(Fresh_cells_filt$Rxrg, Fixed_cells_filt$Rxrg, paired=FALSE)
wilcox.test(Fresh_cells_filt$Arc, Fixed_cells_filt$Arc, paired=FALSE)
wilcox.test(Fresh_cells_filt$Rem2, Fixed_cells_filt$Rem2, paired=FALSE)
wilcox.test(Fresh_cells_filt$Hmbs, Fixed_cells_filt$Hmbs, paired=FALSE)
wilcox.test(Fresh_cells_filt$Sdha, Fixed_cells_filt$Sdha, paired=FALSE)
wilcox.test(Fresh_cells_filt$Gapdh, Fixed_cells_filt$Gapdh, paired=FALSE)
wilcox.test(Fresh_cells_filt$Casp9, Fixed_cells_filt$Casp9, paired=FALSE)
wilcox.test(Fresh_cells_filt$Casp1, Fixed_cells_filt$Casp1, paired=FALSE)




#Means before filtering
#Fresh cells
Fresh_cells_Atf3_mean<-mean(Fresh_cells$Atf3, na.rm=TRUE)
Fresh_cells_Atf3_mean
Fresh_cells_Junb_mean<-mean(Fresh_cells$Junb, na.rm=TRUE)
Fresh_cells_Junb_mean
Fresh_cells_Ier3_mean<-mean(Fresh_cells$Ier3, na.rm=TRUE)
Fresh_cells_Ier3_mean
Fresh_cells_Socs3_mean<-mean(Fresh_cells$Socs3, na.rm=TRUE)
Fresh_cells_Socs3_mean
Fresh_cells_P2ry13_mean<-mean(Fresh_cells$P2ry13, na.rm=TRUE)
Fresh_cells_P2ry13_mean
Fresh_cells_Cd63_mean<-mean(Fresh_cells$Cd63, na.rm=TRUE)
Fresh_cells_Cd63_mean
Fresh_cells_Cldn1_mean<-mean(Fresh_cells$Cldn1, na.rm=TRUE)
Fresh_cells_Cldn1_mean
Fresh_cells_Rxrg_mean<-mean(Fresh_cells$Rxrg, na.rm=TRUE)
Fresh_cells_Rxrg_mean
Fresh_cells_Arc_mean<-mean(Fresh_cells$Arc, na.rm=TRUE)
Fresh_cells_Arc_mean
Fresh_cells_Rem2_mean<-mean(Fresh_cells$Rem2, na.rm=TRUE)
Fresh_cells_Rem2_mean
Fresh_cells_Hmbs_mean<-mean(Fresh_cells$Hmbs, na.rm=TRUE)
Fresh_cells_Hmbs_mean
Fresh_cells_Sdha_mean<-mean(Fresh_cells$Sdha, na.rm=TRUE)
Fresh_cells_Sdha_mean
Fresh_cells_Gapdh_mean<-mean(Fresh_cells$Gapdh, na.rm=TRUE)
Fresh_cells_Gapdh_mean
Fresh_cells_Casp9_mean<-mean(Fresh_cells$Casp9, na.rm=TRUE)
Fresh_cells_Casp9_mean
Fresh_cells_Casp1_mean<-mean(Fresh_cells$Casp1, na.rm=TRUE)
Fresh_cells_Casp1_mean

#Fixed cells
Fixed_cells_Atf3_mean<-mean(Fixed_cells$Atf3, na.rm=TRUE)
Fixed_cells_Atf3_mean
Fixed_cells_Junb_mean<-mean(Fixed_cells$Junb, na.rm=TRUE)
Fixed_cells_Junb_mean
Fixed_cells_Ier3_mean<-mean(Fixed_cells$Ier3, na.rm=TRUE)
Fixed_cells_Ier3_mean
Fixed_cells_Socs3_mean<-mean(Fixed_cells$Socs3, na.rm=TRUE)
Fixed_cells_Socs3_mean
Fixed_cells_P2ry13_mean<-mean(Fixed_cells$P2ry13, na.rm=TRUE)
Fixed_cells_P2ry13_mean
Fixed_cells_Cd63_mean<-mean(Fixed_cells$Cd63, na.rm=TRUE)
Fixed_cells_Cd63_mean
Fixed_cells_Cldn1_mean<-mean(Fixed_cells$Cldn1, na.rm=TRUE)
Fixed_cells_Cldn1_mean
Fixed_cells_Rxrg_mean<-mean(Fixed_cells$Rxrg, na.rm=TRUE)
Fixed_cells_Rxrg_mean
Fixed_cells_Arc_mean<-mean(Fixed_cells$Arc, na.rm=TRUE)
Fixed_cells_Arc_mean
Fixed_cells_Rem2_mean<-mean(Fixed_cells$Rem2, na.rm=TRUE)
Fixed_cells_Rem2_mean
Fixed_cells_Hmbs_mean<-mean(Fixed_cells$Hmbs, na.rm=TRUE)
Fixed_cells_Hmbs_mean
Fixed_cells_Sdha_mean<-mean(Fixed_cells$Sdha, na.rm=TRUE)
Fixed_cells_Sdha_mean
Fixed_cells_Gapdh_mean<-mean(Fixed_cells$Gapdh, na.rm=TRUE)
Fixed_cells_Gapdh_mean
Fixed_cells_Casp9_mean<-mean(Fixed_cells$Casp9, na.rm=TRUE)
Fixed_cells_Casp9_mean
Fixed_cells_Casp1_mean<-mean(Fixed_cells$Casp1, na.rm=TRUE)
Fixed_cells_Casp1_mean


#Means After filtering
#Fresh cells
Fresh_cells_filt_Atf3_mean<-mean(Fresh_cells_filt$Atf3, na.rm=TRUE)
Fresh_cells_filt_Atf3_mean
Fresh_cells_filt_Junb_mean<-mean(Fresh_cells_filt$Junb, na.rm=TRUE)
Fresh_cells_filt_Junb_mean
Fresh_cells_filt_Ier3_mean<-mean(Fresh_cells_filt$Ier3, na.rm=TRUE)
Fresh_cells_filt_Ier3_mean
Fresh_cells_filt_Socs3_mean<-mean(Fresh_cells_filt$Socs3, na.rm=TRUE)
Fresh_cells_filt_Socs3_mean
Fresh_cells_filt_P2ry13_mean<-mean(Fresh_cells_filt$P2ry13, na.rm=TRUE)
Fresh_cells_filt_P2ry13_mean
Fresh_cells_filt_Cd63_mean<-mean(Fresh_cells_filt$Cd63, na.rm=TRUE)
Fresh_cells_filt_Cd63_mean
Fresh_cells_filt_Cldn1_mean<-mean(Fresh_cells_filt$Cldn1, na.rm=TRUE)
Fresh_cells_filt_Cldn1_mean
Fresh_cells_filt_Rxrg_mean<-mean(Fresh_cells_filt$Rxrg, na.rm=TRUE)
Fresh_cells_filt_Rxrg_mean
Fresh_cells_filt_Arc_mean<-mean(Fresh_cells_filt$Arc, na.rm=TRUE)
Fresh_cells_filt_Arc_mean
Fresh_cells_filt_Rem2_mean<-mean(Fresh_cells_filt$Rem2, na.rm=TRUE)
Fresh_cells_filt_Rem2_mean
Fresh_cells_filt_Hmbs_mean<-mean(Fresh_cells_filt$Hmbs, na.rm=TRUE)
Fresh_cells_filt_Hmbs_mean
Fresh_cells_filt_Sdha_mean<-mean(Fresh_cells_filt$Sdha, na.rm=TRUE)
Fresh_cells_filt_Sdha_mean
Fresh_cells_filt_Gapdh_mean<-mean(Fresh_cells_filt$Gapdh, na.rm=TRUE)
Fresh_cells_filt_Gapdh_mean
Fresh_cells_filt_Casp9_mean<-mean(Fresh_cells_filt$Casp9, na.rm=TRUE)
Fresh_cells_filt_Casp9_mean
Fresh_cells_filt_Casp1_mean<-mean(Fresh_cells_filt$Casp1, na.rm=TRUE)
Fresh_cells_filt_Casp1_mean

#Fixed cells_filt
Fixed_cells_filt_Atf3_mean<-mean(Fixed_cells_filt$Atf3, na.rm=TRUE)
Fixed_cells_filt_Atf3_mean
Fixed_cells_filt_Junb_mean<-mean(Fixed_cells_filt$Junb, na.rm=TRUE)
Fixed_cells_filt_Junb_mean
Fixed_cells_filt_Ier3_mean<-mean(Fixed_cells_filt$Ier3, na.rm=TRUE)
Fixed_cells_filt_Ier3_mean
Fixed_cells_filt_Socs3_mean<-mean(Fixed_cells_filt$Socs3, na.rm=TRUE)
Fixed_cells_filt_Socs3_mean
Fixed_cells_filt_P2ry13_mean<-mean(Fixed_cells_filt$P2ry13, na.rm=TRUE)
Fixed_cells_filt_P2ry13_mean
Fixed_cells_filt_Cd63_mean<-mean(Fixed_cells_filt$Cd63, na.rm=TRUE)
Fixed_cells_filt_Cd63_mean
Fixed_cells_filt_Cldn1_mean<-mean(Fixed_cells_filt$Cldn1, na.rm=TRUE)
Fixed_cells_filt_Cldn1_mean
Fixed_cells_filt_Rxrg_mean<-mean(Fixed_cells_filt$Rxrg, na.rm=TRUE)
Fixed_cells_filt_Rxrg_mean
Fixed_cells_filt_Arc_mean<-mean(Fixed_cells_filt$Arc, na.rm=TRUE)
Fixed_cells_filt_Arc_mean
Fixed_cells_filt_Rem2_mean<-mean(Fixed_cells_filt$Rem2, na.rm=TRUE)
Fixed_cells_filt_Rem2_mean
Fixed_cells_filt_Hmbs_mean<-mean(Fixed_cells_filt$Hmbs, na.rm=TRUE)
Fixed_cells_filt_Hmbs_mean
Fixed_cells_filt_Sdha_mean<-mean(Fixed_cells_filt$Sdha, na.rm=TRUE)
Fixed_cells_filt_Sdha_mean
Fixed_cells_filt_Gapdh_mean<-mean(Fixed_cells_filt$Gapdh, na.rm=TRUE)
Fixed_cells_filt_Gapdh_mean
Fixed_cells_filt_Casp9_mean<-mean(Fixed_cells_filt$Casp9, na.rm=TRUE)
Fixed_cells_filt_Casp9_mean
Fixed_cells_filt_Casp1_mean<-mean(Fixed_cells_filt$Casp1, na.rm=TRUE)
Fixed_cells_filt_Casp1_mean

#Summary means
Stress_genes_2<-c("Atf3", "Junb", "Ier3", "Socs3", "P2ry13", "Cd63", "Cldn1", "Rxrg", 
                "Arc", "Rem2", "Hmbs", "Sdha", "Gapdh", "Casp9","Casp1")
Fresh_cells_Before_filt_mean_2<-c(Fresh_cells_Atf3_mean,Fresh_cells_Junb_mean,Fresh_cells_Ier3_mean,Fresh_cells_Socs3_mean,Fresh_cells_P2ry13_mean,
                                 Fresh_cells_Cd63_mean,Fresh_cells_Cldn1_mean,Fresh_cells_Rxrg_mean,Fresh_cells_Arc_mean,Fresh_cells_Rem2_mean,
                                 Fresh_cells_Hmbs_mean,Fresh_cells_Sdha_mean,Fresh_cells_Gapdh_mean,Fresh_cells_Casp9_mean,Fresh_cells_Casp1_mean)

Fixed_cells_Before_filt_mean_2<-c(Fixed_cells_Atf3_mean,Fixed_cells_Junb_mean,Fixed_cells_Ier3_mean,Fixed_cells_Socs3_mean,Fixed_cells_P2ry13_mean,
                                  Fixed_cells_Cd63_mean,Fixed_cells_Cldn1_mean,Fixed_cells_Rxrg_mean,Fixed_cells_Arc_mean,Fixed_cells_Rem2_mean,
                                  Fixed_cells_Hmbs_mean,Fixed_cells_Sdha_mean,Fixed_cells_Gapdh_mean,Fixed_cells_Casp9_mean,Fixed_cells_Casp1_mean)

Fresh_cells_After_filt_mean_2<-c(Fresh_cells_filt_Atf3_mean,Fresh_cells_filt_Junb_mean,Fresh_cells_filt_Ier3_mean,Fresh_cells_filt_Socs3_mean,Fresh_cells_filt_P2ry13_mean,
                                 Fresh_cells_filt_Cd63_mean,Fresh_cells_filt_Cldn1_mean,Fresh_cells_filt_Rxrg_mean,Fresh_cells_filt_Arc_mean,Fresh_cells_filt_Rem2_mean,
                                 Fresh_cells_filt_Hmbs_mean,Fresh_cells_filt_Sdha_mean,Fresh_cells_filt_Gapdh_mean,Fresh_cells_filt_Casp9_mean,Fresh_cells_filt_Casp1_mean)

Fixed_cells_filt_After_filt_mean_2<-c(Fixed_cells_filt_Atf3_mean,Fixed_cells_filt_Junb_mean,Fixed_cells_filt_Ier3_mean,Fixed_cells_filt_Socs3_mean,Fixed_cells_filt_P2ry13_mean,
                                 Fixed_cells_filt_Cd63_mean,Fixed_cells_filt_Cldn1_mean,Fixed_cells_filt_Rxrg_mean,Fixed_cells_filt_Arc_mean,Fixed_cells_filt_Rem2_mean,
                                 Fixed_cells_filt_Hmbs_mean,Fixed_cells_filt_Sdha_mean,Fixed_cells_filt_Gapdh_mean,Fixed_cells_filt_Casp9_mean,Fixed_cells_filt_Casp1_mean)

Stress_genes_summary_mean_2<-data.frame(Stress_genes_2,Fresh_cells_Before_filt_mean_2,Fixed_cells_Before_filt_mean_2,Fresh_cells_After_filt_mean_2,Fixed_cells_filt_After_filt_mean_2)
Stress_genes_summary_mean_2

#standard error of the mean (SEM)

std <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

#stds before filtering
#Fresh cells
Fresh_cells_Atf3_std<-std(Fresh_cells$Atf3, na.rm=TRUE)
Fresh_cells_Atf3_std
Fresh_cells_Junb_std<-std(Fresh_cells$Junb, na.rm=TRUE)
Fresh_cells_Junb_std
Fresh_cells_Ier3_std<-std(Fresh_cells$Ier3, na.rm=TRUE)
Fresh_cells_Ier3_std
Fresh_cells_Socs3_std<-std(Fresh_cells$Socs3, na.rm=TRUE)
Fresh_cells_Socs3_std
Fresh_cells_P2ry13_std<-std(Fresh_cells$P2ry13, na.rm=TRUE)
Fresh_cells_P2ry13_std
Fresh_cells_Cd63_std<-std(Fresh_cells$Cd63, na.rm=TRUE)
Fresh_cells_Cd63_std
Fresh_cells_Cldn1_std<-std(Fresh_cells$Cldn1, na.rm=TRUE)
Fresh_cells_Cldn1_std
Fresh_cells_Rxrg_std<-std(Fresh_cells$Rxrg, na.rm=TRUE)
Fresh_cells_Rxrg_std
Fresh_cells_Arc_std<-std(Fresh_cells$Arc, na.rm=TRUE)
Fresh_cells_Arc_std
Fresh_cells_Rem2_std<-std(Fresh_cells$Rem2, na.rm=TRUE)
Fresh_cells_Rem2_std
Fresh_cells_Hmbs_std<-std(Fresh_cells$Hmbs, na.rm=TRUE)
Fresh_cells_Hmbs_std
Fresh_cells_Sdha_std<-std(Fresh_cells$Sdha, na.rm=TRUE)
Fresh_cells_Sdha_std
Fresh_cells_Gapdh_std<-std(Fresh_cells$Gapdh, na.rm=TRUE)
Fresh_cells_Gapdh_std
Fresh_cells_Casp9_std<-std(Fresh_cells$Casp9, na.rm=TRUE)
Fresh_cells_Casp9_std
Fresh_cells_Casp1_std<-std(Fresh_cells$Casp1, na.rm=TRUE)
Fresh_cells_Casp1_std

#Fixed cells
Fixed_cells_Atf3_std<-std(Fixed_cells$Atf3, na.rm=TRUE)
Fixed_cells_Atf3_std
Fixed_cells_Junb_std<-std(Fixed_cells$Junb, na.rm=TRUE)
Fixed_cells_Junb_std
Fixed_cells_Ier3_std<-std(Fixed_cells$Ier3, na.rm=TRUE)
Fixed_cells_Ier3_std
Fixed_cells_Socs3_std<-std(Fixed_cells$Socs3, na.rm=TRUE)
Fixed_cells_Socs3_std
Fixed_cells_P2ry13_std<-std(Fixed_cells$P2ry13, na.rm=TRUE)
Fixed_cells_P2ry13_std
Fixed_cells_Cd63_std<-std(Fixed_cells$Cd63, na.rm=TRUE)
Fixed_cells_Cd63_std
Fixed_cells_Cldn1_std<-std(Fixed_cells$Cldn1, na.rm=TRUE)
Fixed_cells_Cldn1_std
Fixed_cells_Rxrg_std<-std(Fixed_cells$Rxrg, na.rm=TRUE)
Fixed_cells_Rxrg_std
Fixed_cells_Arc_std<-std(Fixed_cells$Arc, na.rm=TRUE)
Fixed_cells_Arc_std
Fixed_cells_Rem2_std<-std(Fixed_cells$Rem2, na.rm=TRUE)
Fixed_cells_Rem2_std
Fixed_cells_Hmbs_std<-std(Fixed_cells$Hmbs, na.rm=TRUE)
Fixed_cells_Hmbs_std
Fixed_cells_Sdha_std<-std(Fixed_cells$Sdha, na.rm=TRUE)
Fixed_cells_Sdha_std
Fixed_cells_Gapdh_std<-std(Fixed_cells$Gapdh, na.rm=TRUE)
Fixed_cells_Gapdh_std
Fixed_cells_Casp9_std<-std(Fixed_cells$Casp9, na.rm=TRUE)
Fixed_cells_Casp9_std
Fixed_cells_Casp1_std<-std(Fixed_cells$Casp1, na.rm=TRUE)
Fixed_cells_Casp1_std

#stds After filtering
#Fresh cells
Fresh_cells_filt_Atf3_std<-std(Fresh_cells_filt$Atf3, na.rm=TRUE)
Fresh_cells_filt_Atf3_std
Fresh_cells_filt_Junb_std<-std(Fresh_cells_filt$Junb, na.rm=TRUE)
Fresh_cells_filt_Junb_std
Fresh_cells_filt_Ier3_std<-std(Fresh_cells_filt$Ier3, na.rm=TRUE)
Fresh_cells_filt_Ier3_std
Fresh_cells_filt_Socs3_std<-std(Fresh_cells_filt$Socs3, na.rm=TRUE)
Fresh_cells_filt_Socs3_std
Fresh_cells_filt_P2ry13_std<-std(Fresh_cells_filt$P2ry13, na.rm=TRUE)
Fresh_cells_filt_P2ry13_std
Fresh_cells_filt_Cd63_std<-std(Fresh_cells_filt$Cd63, na.rm=TRUE)
Fresh_cells_filt_Cd63_std
Fresh_cells_filt_Cldn1_std<-std(Fresh_cells_filt$Cldn1, na.rm=TRUE)
Fresh_cells_filt_Cldn1_std
Fresh_cells_filt_Rxrg_std<-std(Fresh_cells_filt$Rxrg, na.rm=TRUE)
Fresh_cells_filt_Rxrg_std
Fresh_cells_filt_Arc_std<-std(Fresh_cells_filt$Arc, na.rm=TRUE)
Fresh_cells_filt_Arc_std
Fresh_cells_filt_Rem2_std<-std(Fresh_cells_filt$Rem2, na.rm=TRUE)
Fresh_cells_filt_Rem2_std
Fresh_cells_filt_Hmbs_std<-std(Fresh_cells_filt$Hmbs, na.rm=TRUE)
Fresh_cells_filt_Hmbs_std
Fresh_cells_filt_Sdha_std<-std(Fresh_cells_filt$Sdha, na.rm=TRUE)
Fresh_cells_filt_Sdha_std
Fresh_cells_filt_Gapdh_std<-std(Fresh_cells_filt$Gapdh, na.rm=TRUE)
Fresh_cells_filt_Gapdh_std
Fresh_cells_filt_Casp9_std<-std(Fresh_cells_filt$Casp9, na.rm=TRUE)
Fresh_cells_filt_Casp9_std
Fresh_cells_filt_Casp1_std<-std(Fresh_cells_filt$Casp1, na.rm=TRUE)
Fresh_cells_filt_Casp1_std

#Fixed cells_filt
Fixed_cells_filt_Atf3_std<-std(Fixed_cells_filt$Atf3, na.rm=TRUE)
Fixed_cells_filt_Atf3_std
Fixed_cells_filt_Junb_std<-std(Fixed_cells_filt$Junb, na.rm=TRUE)
Fixed_cells_filt_Junb_std
Fixed_cells_filt_Ier3_std<-std(Fixed_cells_filt$Ier3, na.rm=TRUE)
Fixed_cells_filt_Ier3_std
Fixed_cells_filt_Socs3_std<-std(Fixed_cells_filt$Socs3, na.rm=TRUE)
Fixed_cells_filt_Socs3_std
Fixed_cells_filt_P2ry13_std<-std(Fixed_cells_filt$P2ry13, na.rm=TRUE)
Fixed_cells_filt_P2ry13_std
Fixed_cells_filt_Cd63_std<-std(Fixed_cells_filt$Cd63, na.rm=TRUE)
Fixed_cells_filt_Cd63_std
Fixed_cells_filt_Cldn1_std<-std(Fixed_cells_filt$Cldn1, na.rm=TRUE)
Fixed_cells_filt_Cldn1_std
Fixed_cells_filt_Rxrg_std<-std(Fixed_cells_filt$Rxrg, na.rm=TRUE)
Fixed_cells_filt_Rxrg_std
Fixed_cells_filt_Arc_std<-std(Fixed_cells_filt$Arc, na.rm=TRUE)
Fixed_cells_filt_Arc_std
Fixed_cells_filt_Rem2_std<-std(Fixed_cells_filt$Rem2, na.rm=TRUE)
Fixed_cells_filt_Rem2_std
Fixed_cells_filt_Hmbs_std<-std(Fixed_cells_filt$Hmbs, na.rm=TRUE)
Fixed_cells_filt_Hmbs_std
Fixed_cells_filt_Sdha_std<-std(Fixed_cells_filt$Sdha, na.rm=TRUE)
Fixed_cells_filt_Sdha_std
Fixed_cells_filt_Gapdh_std<-std(Fixed_cells_filt$Gapdh, na.rm=TRUE)
Fixed_cells_filt_Gapdh_std
Fixed_cells_filt_Casp9_std<-std(Fixed_cells_filt$Casp9, na.rm=TRUE)
Fixed_cells_filt_Casp9_std
Fixed_cells_filt_Casp1_std<-std(Fixed_cells_filt$Casp1, na.rm=TRUE)
Fixed_cells_filt_Casp1_std


#Summary stds
Stress_genes_2<-c("Atf3", "Junb", "Ier3", "Socs3", "P2ry13", "Cd63", "Cldn1", "Rxrg", 
                  "Arc", "Rem2", "Hmbs", "Sdha", "Gapdh", "Casp9","Casp1")
Fresh_cells_Before_filt_std_2<-c(Fresh_cells_Atf3_std,Fresh_cells_Junb_std,Fresh_cells_Ier3_std,Fresh_cells_Socs3_std,Fresh_cells_P2ry13_std,
                                  Fresh_cells_Cd63_std,Fresh_cells_Cldn1_std,Fresh_cells_Rxrg_std,Fresh_cells_Arc_std,Fresh_cells_Rem2_std,
                                  Fresh_cells_Hmbs_std,Fresh_cells_Sdha_std,Fresh_cells_Gapdh_std,Fresh_cells_Casp9_std,Fresh_cells_Casp1_std)
Fixed_cells_Before_filt_std_2<-c(Fixed_cells_Atf3_std,Fixed_cells_Junb_std,Fixed_cells_Ier3_std,Fixed_cells_Socs3_std,Fixed_cells_P2ry13_std,
                                  Fixed_cells_Cd63_std,Fixed_cells_Cldn1_std,Fixed_cells_Rxrg_std,Fixed_cells_Arc_std,Fixed_cells_Rem2_std,
                                  Fixed_cells_Hmbs_std,Fixed_cells_Sdha_std,Fixed_cells_Gapdh_std,Fixed_cells_Casp9_std,Fixed_cells_Casp1_std)
Fresh_cells_After_filt_std_2<-c(Fresh_cells_filt_Atf3_std,Fresh_cells_filt_Junb_std,Fresh_cells_filt_Ier3_std,Fresh_cells_filt_Socs3_std,Fresh_cells_filt_P2ry13_std,
                                 Fresh_cells_filt_Cd63_std,Fresh_cells_filt_Cldn1_std,Fresh_cells_filt_Rxrg_std,Fresh_cells_filt_Arc_std,Fresh_cells_filt_Rem2_std,
                                 Fresh_cells_filt_Hmbs_std,Fresh_cells_filt_Sdha_std,Fresh_cells_filt_Gapdh_std,Fresh_cells_filt_Casp9_std,Fresh_cells_filt_Casp1_std)
Fixed_cells_filt_After_filt_std_2<-c(Fixed_cells_filt_Atf3_std,Fixed_cells_filt_Junb_std,Fixed_cells_filt_Ier3_std,Fixed_cells_filt_Socs3_std,Fixed_cells_filt_P2ry13_std,
                                      Fixed_cells_filt_Cd63_std,Fixed_cells_filt_Cldn1_std,Fixed_cells_filt_Rxrg_std,Fixed_cells_filt_Arc_std,Fixed_cells_filt_Rem2_std,
                                      Fixed_cells_filt_Hmbs_std,Fixed_cells_filt_Sdha_std,Fixed_cells_filt_Gapdh_std,Fixed_cells_filt_Casp9_std,Fixed_cells_filt_Casp1_std)
Stress_genes_summary_std_2<-data.frame(Stress_genes_2,Fresh_cells_Before_filt_std_2,Fixed_cells_Before_filt_std_2,Fresh_cells_After_filt_std_2,Fixed_cells_filt_After_filt_std_2)
Stress_genes_summary_std_2



time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Stress_genes_summary_mean_2, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary_2.xlsx"))
           , sheetName = "Stress_genes_summary_mean_2", col.names = TRUE, row.names = TRUE)
write.xlsx(Stress_genes_summary_std_2, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_genes_summary_2.xlsx"))
           , sheetName = "Stress_genes_summary_std_2", col.names = TRUE, row.names = TRUE, append = TRUE)


#Visualization stress genes
DotPlot(object= pilot_filt, features = c("Atf3", "Junb", "Ier3", "Socs3", "P2ry13", "Cd63", "Cldn1", "Rxrg", 
                                         "Arc", "Rem2","Casp9","Casp1","Egr1", "Bax"), dot.scale = 10,scale=FALSE)+  RotatedAxis()
DoHeatmap(object= pilot_filt, slot="counts", features = c("Atf3", "Junb", "Ier3", "Socs3", "P2ry13", "Cd63", "Cldn1", "Rxrg", 
                                           "Arc", "Rem2","Casp9","Casp1","Egr1", "Bax"))

VlnPlot(object= pilot.integrated, features = c("Atf3", "Junb", "Ier3", "Socs3", "P2ry13", "Cd63", "Cldn1", "Rxrg", 
                                         "Arc", "Rem2","Casp9","Casp1","Egr1", "Bax"))
#the boxplots would appear as informative as the violin plots... maybe it isn't the best way of representing at the end



# Wilcoxon Rank Sum test done with Findmarkers function from all the stress genes

Stress_genes_1<-c("Atf3", "Junb", "Ier3", "Socs3", "P2ry13", "Cd63", "Cldn1", "Rxrg", 
                  "Arc", "Rem2","Casp9","Casp1","Egr1", "Bax", "Hsp90aa1","Hsp90ab1", "Hspa1a", "Hspa1b", "Hspa8",
                  "Hspb1", "Hspe1", "Hsph1")
Stress_genes_DE<-FindMarkers(pilot_filt,ident.1 ="Fresh_cells",ident.2 ="Fixed_cells",assay="RNA", slot="counts", features=Stress_genes_1)
Stress_genes_DE

dittoDotPlot(pilot_filt, Stress_genes_1,group.by="Condition",assay="RNA", slot="counts", scale=FALSE, y.reorder=c(2,1),
             max.color =  "#D73027", min.color= "#D0D0D1")

#Seurat object with only the stress genes
#I do a reduced dataset with only the stress genes to do the DE analysis 
#and obtain which ones are differential expressed in fresh
pilot.filt.df<-data.frame(pilot_filt@assays[["RNA"]]@counts)
Pilot.stress<-pilot.filt.df[Stress_genes_1,]

#Export dataframe to test on Sigma Plot and show to Dr Lehmann

time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Pilot.stress, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Pilot.stress.xlsx"))
           , sheetName = "Pilot.stress", col.names = TRUE, row.names = TRUE, append = FALSE)

#Run the DE analysis

Pilot_stress<-CreateSeuratObject(counts = Pilot.stress, project = "Pilot analysis", min.cells = 0, min.features = 0)

current.ids2 <- c("JEN.AU.s004","JEN.AU.s002")
new.ids2 <- c("Fresh_cells","Fixed_cells")
Pilot_stress@active.ident <- plyr::mapvalues(x = Pilot_stress@active.ident, from = current.ids2, to = new.ids2)
Pilot_stress$Condition <- Idents(Pilot_stress)

Pilot_stress$Condition <- factor(x=Pilot_stress$Condition, levels=c("Fresh_cells","Fixed_cells"))




Stress_DE<-FindMarkers(Pilot_stress,ident.1 ="Fresh_cells",ident.2 ="Fixed_cells",assay="RNA", slot="counts")
Stress_DE
Stress_DE$BH_p_val_adj<-p.adjust(Stress_DE$p_val,method="BH")
Stress_DE

dittoDotPlot(Pilot_stress, Stress_genes_1,group.by="Condition",assay="RNA", slot="counts", scale=FALSE, y.reorder=c(2,1),
             max.color =  "#D73027", min.color= "#D0D0D1")






#To obtain the full list of genes from FindMarkers
options(scipen = 999)
Stress_DE_all<-FindMarkers(Pilot_stress,ident.1 ="Fresh_cells",ident.2 ="Fixed_cells",assay="RNA", slot="counts",  
                         min.cells.group = 1, 
                         min.cells.feature = 1,
                         min.pct = 0,
                         logfc.threshold = 0,
                         only.pos = FALSE)
Stress_DE_all
Stress_DE_all$BH_p_val_adj<-p.adjust(Stress_DE_all$p_val,method="BH")
Stress_DE_all


time<- paste(Sys.time())
time.name<-gsub(":","..",time)
write.xlsx(Stress_DE_all, noquote(paste("F:/Marta/Getting started with scRNAseq analysis/Pilot results/Pilot analysis + results/3 Check the quality control parameters from our cells/", paste(time.name), "Stress_DE.xlsx"))
           , sheetName = "Stress_DE", col.names = TRUE, row.names = TRUE, append = FALSE)
