# mESCs-scRNAseq

This repository contains R code and useful functions for analysing mouse embryonic stem cell (mESC) single cell RNAseq data from 
Kolodziejczyk AA et al. 2015. Single cell RNA-sequencing of pluripotent states unlocks modular transcriptional variation. Cell Stem Cell.
17:471–485:  [](http://www.sciencedirect.com/science/article/pii/S193459091500418X).

The data has been processed by Luke Zappia and can be found in `/group/bioi1/shared/public_data/Kolodziejczyk-mESCs-scRNAseq/`. 
The counts table is at `/group/bioi1/shared/public_data/Kolodziejczyk-mESCs-scRNAseq/counts.counts.txt` and the annotation information is
at `/group/bioi1/shared/public_data/Kolodziejczyk-mESCs-scRNAseq/metadata/E-MTAB-2600.sdrf.txt`.

## Directory structure

* **R** - Resuable R code (functions etc.)
* **analysis** - RMarkdown analysis files

## Analysis

Current analysis files include:

* `Filtering-mESCs-01.Rmd` - Quality control and filtering of genes and cells
* `scater_qc.Rmd` - Scater QC plots of the full data set
* `scater_qc_filtered.Rmd` Scater QC plots of the filtered dataset

## Code

Current code files include:

* `load_SCESet.R` - Functions for loading data into and SCESet object
* `filter_SCESet.R` - Functions for filtering an SCESet object
