library(dplyr)
library(ggplot2)
library(ggthemes)
library(magrittr)
library(ggpubr)
library(biomaRt)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(dendextend)
library(tibble)

setwd("C:/Users/Erkut Ilaslan/Desktop/data_integration/NANOS1-3/RNA-Seq OE")

#http://venus.ifca.unican.es/Rintro/dataStruct.html nice basics explained


#testing heatmap visuzalization on normalized counts imported from galaxy deseq2
#https://www.biostars.org/p/317349/#317379 perfect example
#another nice tuto https://sebastianraschka.com/Articles/heatmaps_in_r.html

#data import

N1_DGE_normcont <- read.table(file = 'N1_DESeq2_output_normalized_counts.tabular', sep = '\t', header = TRUE)
N3_DGE_normcont <- read.table(file = 'N3_DESeq2_output_normalized_counts.tabular', sep = '\t', header = TRUE)

#chaning colnames

colnames(N1_DGE_normcont) = c('geneID',
                             'pCMV6 1',
                             'pCMV6 2',
                             'pCMV6 3',
                             'NANOS1 1',
                             'NANOS1 2',
                             'NANOS1 3')

colnames(N3_DGE_normcont) = c('geneID',
                              'pCMV6 1',
                              'pCMV6 2',
                              'pCMV6 3',
                              'NANOS3 1',
                              'NANOS3 2',
                              'NANOS3 3')

#this could be also done by merge fucntion I guess it is also safer
#because i will not need to manually check whether the order of rows are the same so
#annoatation will be correct

NANOS_DGE_normcount <- mutate(N1_DGE_normcont, 
                              N3_DGE_normcont$`NANOS3 1`, 
                              N3_DGE_normcont$`NANOS3 2`, 
                              N3_DGE_normcont$`NANOS3 3`)

colnames(NANOS_DGE_normcount)[8:10] = c('NANOS3 1', 'NANOS3 2', 'NANOS3 3') 

head(NANOS_DGE_normcount)

#geneIDs to HGNC
#---------------------------------------------------------------------------------

#http://www.ensembl.org/Help/ArchiveList
#latest Ensembl GRCh38.p5 assembly is version 91, you can access it via the direct:
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# to take a snapshot and get reproducible results we can access the archive like this:
# use version 85 from an archive
# ensembl85=useMart("ENSEMBL_MART_ENSEMBL", host="jul2016.archive.ensembl.org/biomart/martservice/", dataset="mmusculus_gene_ensembl")

# ??biomaRt
# see what's available:
datasets <- listDatasets(ensembl)
attributes <- listAttributes(ensembl)

# we will map all transcripts from entrez to UniProt
# and then create a subset for significant hits only
# we will match both of those against our protein dataset


# prepare query identifiers from the transcripts dataset
gene_list <- N1_DGE_normcont$geneID
#?getBM

ID_to_HGNC <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                filters = c(filters = "entrezgene_id"),
                values = gene_list,
                uniqueRows = TRUE,
                useCache = TRUE,
                mart = ensembl)

#removing duplicated rows
#?duplicated
ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]

# replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.

ID_to_HGNC[ID_to_HGNC == ""] = NA
ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC),]

#add HGNC to dataframe by joining them

NANOS_DGE_normcount_HGNC <-left_join(NANOS_DGE_normcount, ID_to_HGNC, 
                                     by = c("geneID" = 'entrezgene_id'),
                                     copy = FALSE)

#NA omit otherwise we can't move HGNC to rownames

NANOS_DGE_normcount_HGNC <- na.omit(NANOS_DGE_normcount_HGNC)

#removing row names

NANOS_DGE_normcount_HGNC <- remove_rownames(NANOS_DGE_normcount_HGNC)

#moving HGNC to rownames

NANOS_DGE_normcount_HGNC <- column_to_rownames(NANOS_DGE_normcount_HGNC, var = 'hgnc_symbol')

#removing geneID column

NANOS_DGE_normcount_HGNC <- NANOS_DGE_normcount_HGNC[-1]
head(NANOS_DGE_normcount_HGNC)

#z-score calculation

NANOS_DGE_normcount_HGNC_zscores <- t(scale(t(NANOS_DGE_normcount_HGNC)))

#selecing of genes to visualize e.g. infertility vectors

NANOS1_male_infertility <- read.table('NANOS1_infertility_genes.txt', header = FALSE, sep = ' ')
NANOS3_male_infertility <- read.table('NANOS3_infertility_genes.txt', header = FALSE, sep = ' ')

NANOS1_bruggerman_germcell <- read.table('NANOS1_bruggerman_genes.txt', header = FALSE, sep = ' ')

#gene DLEU1 is excluded because it is not present in our dataframe
NANOS3_bruggerman_germcell <- read.table('NANOS3_bruggerman_genes.txt', header = FALSE, sep = ' ')

N1_male_infertility_genes <- NANOS1_male_infertility$V1
N3_male_infertility_genes <- NANOS3_male_infertility$V1

NANOS1_bruggerman <- NANOS1_bruggerman_germcell$V1 
NANOS3_bruggerman <- NANOS3_bruggerman_germcell$V1

#need to take the columns with the genes i want to visualize

NANOS1_heatmap_infertility <- NANOS_DGE_normcount_HGNC_zscores[as.character(N1_male_infertility_genes), 
                                                               c(-7, -8, -9) ]

NANOS3_heatmap_infertility <- NANOS_DGE_normcount_HGNC_zscores[as.character(N3_male_infertility_genes), 
                                                               c(-4, -5, -6)]

NANOS1_heatmap_bruggerman <- NANOS_DGE_normcount_HGNC_zscores[as.character(NANOS1_bruggerman), 
                                                               c(-7, -8, -9) ]

NANOS3_heatmap_bruggerman <- NANOS_DGE_normcount_HGNC_zscores[as.character(NANOS3_bruggerman), 
                                                               c(-4, -5, -6)]

#visualization
#https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/ComplexHeatmap/inst/doc/ComplexHeatmap.html#toc_0
#adamlar efsane yapmis

#?Heatmap

Heatmap(
  NANOS1_heatmap_infertility,
  column_title = 'NANOS1: Male Infertility',
  heatmap_legend_param = list(title = 'z-score', 
                              legend_height = unit(5, "cm"),
                              grid_width = unit(0.75, "cm"),
                              labels_gp = gpar(fontsize = 12),
                              title_gp = gpar(fontsize = 14)),
  column_title_gp = gpar(fontsize = 18),
  #col = greenred(75), 
  #show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 16),
  #show_column_names = TRUE
  #cluster_rows = FALSE,
  cluster_columns = FALSE,
  #show_column_dend = TRUE,
  show_row_dend = FALSE,
  #row_dend_reorder = TRUE,
  #column_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = 'pearson',
  #clustering_method_columns = "ward.D2",
  width = unit(150, "mm"),
  height = unit(130, 'mm')
)

Heatmap(
  NANOS3_heatmap_infertility,
  column_title = 'NANOS3: Male Infertility',
  column_title_gp = gpar(fontsize = 18),
  heatmap_legend_param = list(title = 'z-score', 
                              legend_height = unit(5, "cm"),
                              grid_width = unit(0.75, "cm"),
                              labels_gp = gpar(fontsize = 12),
                              title_gp = gpar(fontsize = 14)),
  #col = greenred(75), 
  #show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 16),
  #show_column_names = TRUE
  #cluster_rows = FALSE,
  cluster_columns = FALSE,
  #show_column_dend = TRUE,
  show_row_dend = FALSE,
  #row_dend_reorder = TRUE,
  #column_dend_reorder = TRUE,
  clustering_distance_rows = 'pearson',
  clustering_method_rows = "ward.D2",
  #clustering_method_columns = "ward.D2",
  width = unit(150, "mm"),
  height = unit(240, 'mm')
)

Heatmap(
  NANOS1_heatmap_bruggerman,
  column_title = 'NANOS1: Germ cell development genes',
  heatmap_legend_param = list(title = 'z-score', 
                              legend_height = unit(5, "cm"),
                              grid_width = unit(0.75, "cm"),
                              labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 16)),
  column_title_gp = gpar(fontsize = 20),
  #col = greenred(75), 
  #show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 18),
  #show_column_names = TRUE
  #cluster_rows = FALSE,
  cluster_columns = FALSE,
  #show_column_dend = TRUE,
  show_row_dend = FALSE,
  #row_dend_reorder = TRUE,
  #column_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = 'pearson',
  #clustering_method_columns = "ward.D2",
  width = unit(150, "mm"),
  height = unit(190, 'mm')
)

Heatmap(
  NANOS3_heatmap_bruggerman,
  column_title = 'NANOS3: Germ cell development genes',
  heatmap_legend_param = list(title = 'z-score', 
                              legend_height = unit(5, "cm"),
                              grid_width = unit(0.75, "cm"),
                              labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 16)),
  column_title_gp = gpar(fontsize = 20),
  #col = greenred(75), 
  #show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 18),
  #show_column_names = TRUE
  #cluster_rows = FALSE,
  cluster_columns = FALSE,
  #show_column_dend = TRUE,
  show_row_dend = FALSE,
  #row_dend_reorder = TRUE,
  #column_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = 'pearson',
  #clustering_method_columns = "ward.D2",
  width = unit(150, "mm"),
  height = unit(360, 'mm')
)

#exporting dimensions 750 x 1200 (for NANOS3_bruggerman 750 x 1800)

