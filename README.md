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
```r
maplot_dge(data, #your data
           FDR = 0.05, #set FDR treshold.
           FC = 0.5, #set fold change threshold.
           TOP = 10, #top significant genes to be labeled on the plot.
           type = "deseq2", #set data type. deseq2 or edger.
           header = "biovizR is amazing!") #set title for the plot.
           
```
<p align="center">
  <img src="https://github.com/erkutilaslan/biovizR/blob/devel/test.jpg" width="450" height="400"></div>
</p>


- Analysis and visualization of qPCR
<p align="center">
  <img src="link_here" width="760" height="412"></div>
</p>

- Barplot visualization of GO Analysis
<p align="center">
  <img src="link_here" width="760" height="412"></div>
</p>

- Heatmap visualization
<p align="center">
  <img src="link_here" width="760" height="412"></div>
</p>

- Violinplot visualization
<p align="center">
  <img src="link_here" width="760" height="412"></div>
</p>

- Dotplot visualization of IP-MS
<p align="center">
  <img src="link_here" width="760" height="412"></div>
</p>

# Examples



