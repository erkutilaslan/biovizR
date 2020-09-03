#how to select subset of genes with names from the whole RNA-Seq. I need to try this later
#for now I will just use cytoscape

library(ggplot2)
library(ggthemes)
library(dplyr)
library(magrittr)
library(ggpubr)
library(biomaRt)
library(tibble)

#setwd('~/data_integration/NANOS1-3/RNA-Seq OE/')
setwd("C:/Users/Erkut Ilaslan/Desktop/data_integration/NANOS1-3/RNA-Seq OE/Workspace")

#N3 downregulated differentiation genes visualization for Kamila NANOS3 grant  

#Importing data 
DGE_normcont <- read.table(file = 'N3_DESeq2_output_normalized_counts.tabular', 
                              sep = '\t', header = TRUE) 

#changing column names colnames
colnames(DGE_normcont) = c('geneID', 'pCMV6 1', 'pCMV6 2', 'pCMV6 3', 'NANOS3 1', 'NANOS3 2', 'NANOS3 3')  

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
gene_list <- DGE_normcont$geneID
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

DGE_normcount_HGNC <-left_join(DGE_normcont, ID_to_HGNC, 
                                     by = c("geneID" = 'entrezgene_id'),
                                     copy = FALSE)

#NA omit otherwise we can't move HGNC to rownames

DGE_normcount_HGNC <- na.omit(DGE_normcount_HGNC)

#removing row names

DGE_normcount_HGNC <- remove_rownames(DGE_normcount_HGNC)

#moving HGNC to rownames

DGE_normcount_HGNC <- column_to_rownames(DGE_normcount_HGNC, var = 'hgnc_symbol')

#removing geneID column

DGE_normcount_HGNC <- DGE_normcount_HGNC[-1]
head(DGE_normcount_HGNC)

#table with differentiation genes  
selected_genes <- c('TFAP2A', 'TFAP2C', 'SOX17', 'POU5F1', 'PRDM14', 
                                  'SALL4', 'NANOG', 'FGF4', 'KLF5', 'KLF4', 'DNMT3L', 
                                  'TBX3', 'GATA6')       

#subsetting differentiation genes
selected_genes_boxplot <- DGE_normcount_HGNC[as.character(selected_genes), ]

#visualization
#i dont know how to visualize this to be honest seems like a stupid idea.

?ggboxplot


