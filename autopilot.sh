#!/bin/bash

# This is the main script to run the entire single-cell RNA-seq data analysis pipeline.
# Please put raw sequencing reads *.fastq.gz files in their corresponding "data/raw" folders before running this script.
# Please make sure you have installed all the required software and packages as specified in the README file.

conda activate main

# run sequence alignment and generate gene-cell matrix
cd cellranger
bash cellranger.sh # run Cell Ranger pipeline
python metrics.py # plot quality control metrics

# preprocessing, quality control, normalization, and clustering
cd ../qc
python qc.py

# plot hemoglobin and calprotectin expression levels before and after filtering
cd ../contamination
python contamination.py # obtain hemoglobin and calprotectin expression levels
python plot.py # calculate metrics and plot expression levels
python cellxgene.py # visualize expression levels from reference atlas

# cell type annotation
cd ../annotation
python signatures.py # output cell type signatures obtained from quality control step
python markers.py # annotate cell types based on known markers and output markers of specific cell types
python plot.py # plot cell type annotation results

# differential expression analysis
cd ../deg
python deg.py # perform differential expression analysis by gravity
python plot.py # plot differential expression results
python summary.py # plot differential expression profiles for selected genes

# accelerated ageing test
cd ../ageing
python ageing.py # perform differential expression analysis by age
python plot.py # plot ageing signatures and correlations

# gene set enrichment analysis
conda activate gsea # GSEA requires a different environment due to conflicts with scanpy and python 3.13
cd ../gsea
python gsea.py # perform gene set enrichment analysis
python plot.py # plot gene set enrichment results
python summary.py # summarize gene set enrichment results with respect to selected pathways
conda deactivate

# causal structure learning
cd ../causal
python discretize.py # discretize gene expression data
Rscript causal.r # perform causal structure learning

# human-mouse comparison analysis
conda activate comparison # Comparison analysis requires specific Python and R dependencies
cd ../comparison
python convert_orthologs.py # map mouse genes to human orthologs
Rscript fisher_test.R # perform Fisher's exact test
Rscript overlap_heatmap.R # plot concordant genes heatmap
conda deactivate

conda deactivate
