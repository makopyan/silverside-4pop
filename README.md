
# Genetic differentiation is constrained to chromosomal inversions and putative centromeres in locally adapted populations with higher gene flow

This repository contains code and documentation for analyzing population genetic structure and genomic features in four Atlantic silverside (*Menidia menidia*) populations along the North American east coast.


## Data 

[Raw reads](http://ncbi.nlm.nih.gov/bioproject/PRJNA376564/)
[Reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_965154125.1/)

## Repository Contents

- `plot.Rmd`: R markdown for generating visualizations of population genetic statistics and genomic features
- `data.Rmd`: R markdown for data processing, SNP calling, and calculation of population genetic statistics
- `/scripts`: Shell scripts for ANGSD analyses, recombination estimation, and feature detection
- `/data_files`: Input data files used by the R markdown code


## Analysis Overview

This project analyzes whole-genome resequencing data from Atlantic silverside populations to examine:
1. Population structure and differentiation (PCA, FST, DXY)
2. Genetic diversity (θ and π)
3. Linkage disequilibrium patterns
4. Inversions and other structural features
5. Relationship between recombination rate and genetic diversity/differentiation

## Data Structure

The analysis follows these general steps:
1. SNP calling with ANGSD
2. Calculating population genetic statistics in 50kb windows
3. Identifying genomic features (inversions, centromeres, telomeres)
4. Integrating genetic and genomic data for visualization and statistical analysis
