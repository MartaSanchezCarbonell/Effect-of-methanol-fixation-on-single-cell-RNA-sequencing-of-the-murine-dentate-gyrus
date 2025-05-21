# 🧠 Effect of Methanol Fixation on Single-Cell RNA-seq of the Murine Dentate Gyrus

This repository contains the full analysis pipeline for the study:

> **Sánchez-Carbonell, M., Jiménez Peinado, P., Bayer-Kaufmann, C., Hennings, J. C., Hofmann, Y., Schmidt, S., Witte, O. W., & Urbach, A. (2023).**  
> *Effect of methanol fixation on single-cell RNA sequencing of the murine dentate gyrus.*  
> Frontiers in Molecular Neuroscience, 16, 1223798.  
> [https://doi.org/10.3389/fnmol.2023.1223798](https://doi.org/10.3389/fnmol.2023.1223798)

---

## 🧬 Overview

Methanol (MeOH) fixation is a promising method for preserving cells prior to scRNA-seq, but its effects on adult neural tissue are underexplored. This study evaluates how MeOH fixation impacts transcriptomic quality and biological interpretation using SORT-seq on murine dentate gyrus samples.

This repository includes all scripts used to:

- Preprocess scRNA-seq data from fresh and MeOH-fixed cells, obtained by SOrting and Robot-assisted Transcriptome SEQuencing (SORT-seq), a partially robotized version of the CEL-seq2 protocol.
- Perform clustering and annotation
- Evaluate fixation-induced technical and biological effects (e.g., ambient RNA, dropout, stress markers, transcript length bias)

---

## 🗂️ Repository Structure
Pre-processing SORT-seq/ # Quality control, filtering, normalization, clustering, annotation

Downstream analysis SORT-seq/ # Stress genes, ambient RNA, dropouts, transcript length, clustering


---

## 🔧 Preprocessing Summary

Using **R v4.2.0** and **Seurat v4.1.0**, preprocessing includes:

- Calculation of quality metrics:
  - Transcripts/cell, genes/cell, % mitochondrial, % ERCC, transcript/gene ratio
- Filtering criteria:
  - 800–35,000 transcripts
  - more than 500 genes
  - more than 500 ERCC counts
  - Transcript-to-gene ratio >1.2
- Normalization with **SCTransform**
- PCA, UMAP, clustering (`FindClusters()` + K-means `k=10`)
- Marker gene expression for cell-type annotation

---

## 🔍 Downstream Analyses

Located in `Downstream analysis SORT-seq/`:

**📈 Cell Stress Signature**
   - Expression of stress-related genes (e.g. *Atf3, Hsp90ab1, Socs3*)

**🧪 Ambient RNA Estimation**
   - Using `DecontX` to assess RNA contamination

**🧬 Dropout Quantification**  
   *(by Patricia Jiménez Peinado)*
   - Dropouts measured across increasing gene expression thresholds

**📏 Transcript Length Bias**
   - Expression changes binned by exon length from GRCm38 annotation

**🧠 Correlation & Clustering**
   - Pairwise correlation matrices and hierarchical clustering of fresh vs fixed cells

---

## 🛠️ Requirements

- **R ≥ 4.2.0**
- R packages:
  - `Seurat`
  - `DecontX`
  - `GenomicFeatures`
  - `AnnotationDbi`
  - `pheatmap`
  - `ggplot2`, `dplyr`, `data.table` etc.

---

## 🚀 Getting Started
**Clone the repository**

git clone https://github.com/MartaSanchezCarbonell/Effect-of-methanol-fixation-on-single-cell-RNA-sequencing-of-the-murine-dentate-gyrus.git
cd Effect-of-methanol-fixation-on-single-cell-RNA-sequencing-of-the-murine-dentate-gyrus

Run scripts in Pre-processing SORT-seq/ to clean and annotate the dataset.

Proceed with Downstream analysis SORT-seq/ scripts in sequence or explore them individually.

---

## 📚 Citation

If you use this repository or its scripts, please cite:

Sánchez-Carbonell, M., et al. (2023). Effect of methanol fixation on single-cell RNA sequencing of the murine dentate gyrus. Frontiers in Molecular Neuroscience, 16, 1223798. https://doi.org/10.3389/fnmol.2023.1223798

---

## 👩‍🔬 Authors
Marta Sánchez-Carbonell – Analysis pipeline and manuscript

Patricia Jiménez Peinado – Dropout evaluation

With contributions from: C. Bayer-Kaufmann, J. C. Hennings, Y. Hofmann, S. Schmidt, O. W. Witte, A. Urbach

---

## 📬 Contact
For questions or collaborations:

GitHub: @MartaSanchezCarbonell

Email: sanchezcarbonellmarta@gmail.com
