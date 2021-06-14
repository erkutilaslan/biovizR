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
# Data Import

biovizR supports both providing a path to the data and importing the data to R Studio via GUI.


```r
barplot_qpcr("C://Users/Erkut/Desktop/qPCR/siFOXM1_qpcr.csv",
             group1 = "siCTRL",
             group2 = "siFOXM1",
             ref1 = "GAPDH",
             ref2 = "ACTB",
             goi = "FOXM1")
             
barplot_qpcr(siFOXM1,
             group1 = "siCTRL",
             group2 = "siFOXM1",
             ref1 = "GAPDH",
             ref2 = "ACTB",
             goi = "FOXM1")
```
Both of these work.

**Importing data using RStudio GUI:**

A more user friendly way to import your data into RStudio is to use the graphical user interface (GUI) elements.

- Step 1:

The enviroment window in RStudio is where the data will be stored as an object. 
Clicking on the "Import Data" button to open will show options for importing different types of data.
<p align="left">
  <img src="https://github.com/erkutilaslan/biovizR/blob/devel/import1.png" width="1700" height="1100"></div>
</p>

- Step 2:

Let's try to import our RT-qPCR results as an excel data. Clicking on the "From Excel..." button will open a new window. 
Here clicking on the "Browse..." button will let you select your data.
<p align="left">
  <img src="https://github.com/erkutilaslan/biovizR/blob/devel/import2.png" width="1700" height="1100"></div>
</p>

- Step 3:

Selecting your data allows you to preview it. Here you can set the name of the object that will contain your data. By default the name of your excel data will be used.
<p align="left">
  <img src="https://github.com/erkutilaslan/biovizR/blob/devel/import3.png" width="1700" height="1100"></div>
</p>

Once you import your data to R enviroment, you can now apply the functions of biovizR to analyze and visualize your data.

Example:
```r
barplot_qpcr(siFOXM1,
             group1 = "siCTRL",
             group2 = "siFOXM1",
             ref1 = "GAPDH",
             ref2 = "ACTB",
             goi = "FOXM1")
```

# Usage

- **MA-plot visualization of RNA-Seq**

This function enables you to visualize differential gene expression data-sets (for example from RNA-Seq) as an MA-plot.

```r
maplot_dge(data, #your data
           FDR = 0.05, #set FDR treshold.
           FC = 1, #set fold change threshold.
           TOP = 10, #top significant genes to be labeled on the plot.
           header = "MA-plot is amazing!", #title for the MA-plot.
           type = "deseq2", #set data type. deseq2 or edger.
           header = "MA-plot is amazing!") #set title for the plot.      
```

<p align="left">
  <img src="https://github.com/erkutilaslan/biovizR/blob/devel/maplot_dge.png" width="500" height="500"></div>
</p>


- **Analysis and visualization of qPCR**
```r
barplot_qpcr(data, #your data
             group1 = "siCTRL", #name of group 1.
             group2 = "siPUM1", #name of group 2.
             ref = "ACTB", #name of the first reference gene.
             ref2 = "GAPDH", #name of the second reference gene if used.
             goi = "NANOS1", #name of gene of interest
             test = TRUE, #set to FALSE if statistical analysis is not applicable.
             stat = "t-test") #statistical analysis method. t-test or anova can be used.
```


<p align="left">
  <img src="https://github.com/erkutilaslan/biovizR/blob/devel/barplot_qpcr.png" width="450" height="450"></div>
</p>

- **Barplot visualization of GO Analysis**

This function visualize gene ontology (GO) results as a barplot.
<p align="left">
  <img src="https://github.com/erkutilaslan/biovizR/blob/devel/barplot_go.png" width="1000" height="450"></div>
</p>

- **Heatmap visualization**

This function scale expression values with z-score and visualzie them as a heatmap.
<p align="left">
  <img src="link_here" width="760" height="412"></div>
</p>

- **Violinplot visualization**

This function visualize biological data as a violin plot.
<p align="left">
  <img src="link_here" width="760" height="412"></div>
</p>

- **Dotplot visualization of IP-MS**

This function visualize immunoprecipitation followed by mass spectrometry results as a dot-plot. 
<p align="left">
  <img src="link_here" width="760" height="412"></div>
</p>



