# biovizR: An R package to analyze and visualize biological data

biovizR is an R package developed for wet lab scientists to provide easy to use tools for analysis and visualization of biological data with minimal programming background.

This package  relies on already available packages like ggplot2, ggpubr and ComplexHeatmap for visualization. If you are using biovizR for data visualization and have knowledge of R programming, it is highly recommended to use these packages directly for more flexibility and customization.


# Requirements

- R version 3.6 or higher

- devtools package

```r
install.packages("devtools")
```

- tidyverse package
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("tidyverse")
```
- ggpubr package
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggpubr")
```
- ComplexHeatmap package
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
```
- biomaRt package
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")
```
# Installation

```r
devtools::install_github("erkutilaslan/biovizR")
```
# Usage

- MA-plot visualization of RNA-Seq

- Analysis and visualization of qPCR

- Barplot visualization of GO Analysis

- Heatmap visualization

- Violinplot visualization

- Dotplot visualization of IP-MS

# Examples



