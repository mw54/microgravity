# Post-spaceflight analysis reveals enduring reprogramming of the murine immune system by microgravity
 
This repository contains analysis scripts for investigating spaceflight effects on immune cells using single-cell RNA sequencing data from NASA's RRRM-2 mission.
 
## Data Source
 
The analysis uses scRNA-seq data from four datasets (GLDS-402, GLDS-403, GLDS-404, GLDS-405) from NASA's RRRM-2 mission, examining different tissues:
- [**402**](https://osdr.nasa.gov/bio/repo/data/studies/OSD-402): Femur bone marrow
- [**403**](https://osdr.nasa.gov/bio/repo/data/studies/OSD-403): Humerus bone marrow
- [**404**](https://osdr.nasa.gov/bio/repo/data/studies/OSD-404): Peripheral blood mononuclear cells (PBMC)
- [**405**](https://osdr.nasa.gov/bio/repo/data/studies/OSD-405): Spleen
 
The datasets were obtained from NASA Space Biology Open Science Data Repository (OSDR) in August 2025 from https://registry.opendata.aws/nasa-osdr. Each dataset compares samples from mice exposed to microgravity (μG) versus ground control (1G), across young and old age groups.
 
## Repository Structure
```
.
├── autopilot.sh                 # Main pipeline execution script
├── cellranger/                  # Sequence alignment and QC
│   ├── cellranger.sh
│   ├── metrics.py
│   └── constants.py
├── qc/                          # Preprocessing and clustering
│   ├── qc.py
│   └── constants.py
├── contamination/               # Contamination assessment
│   ├── data
│   │   └── cellxgene.csv
│   ├── contamination.py
│   ├── plot.py
│   ├── cellxgene.py
│   └── constants.py
├── annotation/                  # Cell type annotation
│   ├── signatures.py
│   ├── markers.py
│   ├── plot.py
│   └── constants.py
├── deg/                         # Differential expression analysis
│   ├── deg.py
│   └── plot.py
├── gsea/                        # Gene set enrichment analysis
│   ├── gsea.py
│   ├── plot.py
│   └── summary.py
├── causal/                      # Causal structure learning
│   ├── discretize.py
│   └── causal.r
└── comparison/                  # Cross-species comparison
    ├── convert_orthologs.py
    ├── fisher_test.R
    └── overlap_heatmap.R
```
 
## Script Overview
 
### Pipeline Execution
- **autopilot.sh**: Automated execution of the complete analysis pipeline from raw sequencing reads to causal structure learning
 
### Cell Ranger Processing
- **cellranger.sh**: Runs Cell Ranger pipeline for sequence alignment and gene-cell matrix generation
- **metrics.py**: Generates quality control metrics and visualization plots for each sample
 
### Quality Control and Preprocessing
- **qc.py**: Performs cell/gene filtering, doublet removal, normalization, dimensionality reduction (PCA/UMAP), and clustering using Leiden algorithm
 
### Contamination Analysis
- **contamination.py**: Assesses hemoglobin and calprotectin contamination levels across droplet types
- **plot.py**: Visualizes contamination markers before and after filtering
- **cellxgene.py**: Visualizes expression levels of relevant genes from CELLxGENE reference atlas
 
### Cell Type Annotation
- **signatures.py**: Extracts cluster-specific gene signatures from quality-controlled data
- **markers.py**: Annotates cell types based on known marker genes and identifies cell type-specific markers
- **plot.py**: Generates UMAP plots, marker dotplots, and cell type proportion visualizations
 
### Differential Expression Analysis
- **deg.py**: Identifies differentially expressed genes between microgravity and ground control conditions for each cell type and age group
- **plot.py**: Creates volcano plots and dotplots summarizing top DEGs
 
### Gene Set Enrichment Analysis
- **gsea.py**: Performs GSEA using MSigDB canonical pathways to identify enriched biological processes
- **plot.py**: Visualizes enrichment results with barplots and dotplots
- **summary.py**: Summarizes enrichment across selected pathway categories (mitochondria, proliferation, motility, immunity)
 
### Causal Structure Learning
- **discretize.py**: Discretizes gene expression data into categorical levels for causal inference
- **causal.r**: Performs Bayesian network learning using bootstrapped Markov blanket discovery to identify genes causally related to gravity conditions
 
### Cross-Species Comparison
- **convert_orthologs.py**: Maps mouse genes to human orthologs using MyGene.info
- **fisher_test.R**: Performs Fisher's exact test to assess statistical enrichment of orthologous genes between human and mouse datasets
- **overlap_heatmap.R**: Visualizes concordant gene expression changes between human and mouse using heatmaps
 
## Requirements
 
- Cell Ranger (v9.0.1)
- Conda
 
### Environment Setup
 
This pipeline requires three separate conda environments due to package compatibility issues:
 
1. **Main environment** (Python 3.13.5) - for most analyses
2. **GSEA environment** (Python 3.12.11) - for gene set enrichment analysis
3. **Comparison environment** (Python 3.12 + R) - for human-mouse comparison
 
Create the environments:
```bash
# Main environment (cellranger, qc, contamination, annotation, deg, causal)
conda env create -f setup/main.yaml
 
# GSEA environment (gsea only)
conda env create -f setup/gsea.yaml

# Comparison environment (comparison only)
conda env create -f setup/comparison.yaml
```
 
## Usage
 
Place raw FASTQ files in their corresponding `data/raw` folders, then execute:
```bash
bash autopilot.sh
```
 
The pipeline will sequentially process all datasets through alignment, QC, annotation, differential expression, enrichment analysis, and causal inference.
 
## License
 
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
