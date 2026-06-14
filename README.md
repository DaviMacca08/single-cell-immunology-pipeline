End-to-End Single-Cell RNA-seq Analysis Workflow: Seurat, Harmony,
SingleR, Azimuth and Consensus Cell-Type Annotation
================

knitr::opts_chunk\$set(echo = TRUE, warning = FALSE, message = FALSE)

![R](https://img.shields.io/badge/R-4.6.0-blue)
![scRNA–seq](https://img.shields.io/badge/scRNA--seq-Seurat-purple)
![Status](https://img.shields.io/badge/status-active%20development-orange)

------------------------------------------------------------------------

# 💼 Bioinformatics Service Demonstration

This repository presents a complete end-to-end single-cell RNA-seq
(scRNA-seq) analysis workflow implemented in R using the Seurat
ecosystem.

It is designed as a reproducible and modular pipeline, suitable for:

- exploratory analysis of scRNA-seq datasets
- quality control and preprocessing
- batch correction and integration
- unsupervised identification of transcriptional populations
- reference-based and marker-supported cell-type annotation

**Analysis capabilities demonstrated in this project include:**

- Quality control assessment (gene counts, UMI counts, mitochondrial
  content)
- Normalization and highly variable gene identification
- Dimensionality reduction (PCA, UMAP)
- Graph-based clustering
- Batch integration using Harmony
- Automated cell-type annotation using reference-based approaches
- Consensus cell-type assignment integrating multiple annotation
  strategies
- Reproducible visualization and reporting

**📤 Deliverables**

- Modular R analysis pipeline
- Integrated Seurat object
- QC and integration diagnostic plots
- Global clustering analysis
- Consensus cell-type annotation at whole-dataset level

# 📊 Case Study

This workflow is demonstrated using a human intestinal single-cell
RNA-seq dataset:

**Organism:** Homo sapiens (human)

**Tissue/Model:** Intestinal mucosa from Crohn’s disease patients and
healthy controls, including lamina propria and intraepithelial
compartments

**Platform:** 10x Genomics Chromium single-cell RNA-seq

**Public dataset source:** GEO (Gene Expression Omnibus) associated with
the study

**Study DOI:** <https://doi.org/10.1038/s41467-021-22164-6>

------------------------------------------------------------------------

# 🧬 Biological Context

The aim of this analysis is to characterize cellular heterogeneity
within intestinal immune compartments by:

- identifying transcriptionally distinct cell populations
- integrating multiple biological samples
- assigning biologically meaningful cell identities
- generating a validated cellular reference for downstream
  compartment-specific analyses

The current analysis focuses on global dataset characterization before
biological subsetting.

------------------------------------------------------------------------

# 📌 Key Results

The global analysis workflow achieved:

- Identification of transcriptionally distinct clusters using
  graph-based clustering
- Correction of sample-associated variation using Harmony integration
  with sample_id as integration variable
- Visualization of integrated cellular structure using PCA and UMAP
- Reference-based annotation using:
  1)  Azimuth mapping
  2)  SingleR reference annotation
  3)  canonical marker validation

A consensus cell-type annotation strategy was applied by integrating
computational predictions with biological marker evaluation.

The resulting annotated Seurat object represents the global cellular
landscape and serves as the starting point for downstream
compartment-specific analyses (e.g. IEL and LPL subsetting).

------------------------------------------------------------------------

# ⚙️ Analysis Pipeline

1.  Data loading & preprocessing

- Raw count matrix import
- Sample metadata integration
- Initial Seurat object creation

2.  Quality control

Quality filtering based on:

- nFeature_RNA
- nCount_RNA
- mitochondrial RNA percentage (percent.mt)

Diagnostic visualization:

- QC violin plots
- feature-count distributions
- mitochondrial content assessment

3.  Normalization and feature selection

- LogNormalize workflow
- Identification of highly variable genes (HVGs)

4.  Dimensionality reduction

- PCA
- UMAP embedding

5.  Batch integration

- Harmony integration

Batch variable:

- sample_id

Integration quality assessed through:

- PCA/UMAP visualization before and after correction
- sample distribution within integrated embeddings

6.  Unsupervised clustering

Graph-based clustering performed on integrated representations.

Outputs:

- transcriptional clusters
- cluster markers
- global cellular landscape

7.  Cell-type annotation

Consensus annotation workflow:

- Azimuth reference mapping
- SingleR automated annotation
- Canonical marker-based validation
- Consensus biological interpretation

## 8. Visualization

The workflow produces a comprehensive set of visualization outputs for
quality assessment, integration evaluation, clustering interpretation
and cell-type annotation.

Selected representative outputs are included in this repository:

- QC distribution plots
- PCA and UMAP embeddings
- Cluster-level visualization
- Reference-based annotation plots (Azimuth / SingleR)
- Feature plots for canonical markers
- Marker expression visualization

Additional plots generated during intermediate analytical steps are
retained locally but not uploaded to the repository to keep the project
lightweight and focused on key results.

------------------------------------------------------------------------

# Experimental Design

- Multi-sample human intestinal scRNA-seq dataset

- Crohn’s disease vs healthy control samples

- Multiple patient-derived biopsies

- Variable number of cells per sample

The global annotation step precedes compartment-specific analyses.

------------------------------------------------------------------------

# 📁 Project Structure

``` text
scRNAseq_project/
├── scripts/
│   ├── 00_setup/
│   ├── 01_QC_single_sample.R
│   ├── 02_integration_harmony.R
│   ├── 03_global_exploration.R
│   ├── 04_cell_annotation_singler.R
│   ├── 05_cell_annotation_azimuth.R
│   ├── 06_consensus_celltype_annotation
│
├── plots/
│   ├── Azimuth_PBMC/
│   ├── CellAnnotation/
│   ├── Global_Exploration/
│   ├── Integration_HVG/
│   ├── QC_by_sample/
│   ├── SingleR/
│   
├── reports/
├── tables/
│ 
├── README.md
```

# References

- Stuart et al., 2019 – Comprehensive Integration of Single-Cell Data
  (Seurat)
- Korsunsky et al., 2019 – Harmony: Fast, sensitive and accurate
  integration
- Hao et al., 2021–2023 – Seurat v4/v5 framework and multimodal
  integration
- Satija Lab – Seurat documentation: <https://satijalab.org/seurat/>
- Azimuth reference mapping: <https://azimuth.hubmapconsortium.org/>

# 🛠 Tools

- R v4.6.0
- Seurat v5.5.0
- Harmony v2.0.3
- SingleR v2.14.0
- Azimuth v0.5.1
- celldex v1.22.0

------------------------------------------------------------------------

# Note

This repository is structured as a reproducible bioinformatics workflow
rather than a single analysis script.

The pipeline emphasizes:

- modularity
- reproducibility
- transparent annotation strategy
- biological interpretability
- suitability for real-world scRNA-seq analysis projects
