#dependencies

library(ggrepel)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(magrittr)
library(ggpubr)
library(biomaRt)
library(readxl)

update.packages('ggpubr')

#setting working directory
setwd('C:/Users/Erkut Ilaslan/Desktop/data_integration/NANOS1-3/RNA-Seq OE/Workspace/')

# loading data
#these data didnt work for this visualziation

#N3_up <- read.table(file = 'NANOS3_up_0.01_HGNC.csv', 
#                    header = TRUE,
#                    sep = ',',
#                    dec = '.') 
#N3_down <- read.table(file = 'NANOS3_down_0.01_HGNC.csv',
#                      header = TRUE,
#                      sep = ',',
#                      dec = '.')
#N1_up <- read.table(file = 'NANOS1_up_0.01_HGNC.csv',
#                    header = TRUE,
#                    sep = ',',
#                    dec = '.')
#N1_down <- read.table(file = 'NANOS1_down_0.01_HGNC.csv',
#                      header = TRUE,
#                      sep = ',',
#                      dec = '.')
#col_names <- c('geneID',	'BaseMean',	'log2FC',	'StdErr',	'Wald-Stats',	'P-value',	'P-adj')

#the following data is the proper one

#DEG_NANOS1 <- read.table(file = 'NANOS1 OE DE.tabular',
#                         header = FALSE)
#colnames(DEG_NANOS1) = col_names
#DEG_NANOS3 <- read.table(file = 'NANOS3 OE DE.tabular',
#                         header = FALSE)
#colnames(DEG_NANOS3) = col_names

#removing unnecessary columns
#again this step is useless this data was not proper for this visualization
#N1_down$SUID <- NULL 
#N1_down$selected <- NULL
#N1_down$P.value <- NULL
#N1_down$shared.name <- NULL

#N1_up$SUID <- NULL
#N1_up$selected <- NULL
#N1_up$P.value <- NULL
#N1_up$shared.name <- NULL

#N3_down$SUID <- NULL
#N3_down$selected <- NULL
#N3_down$P.value <- NULL
#N3_down$shared.name <- NULL

#N3_up$SUID <- NULL
#N3_up$selected <- NULL
#N3_up$P.value <- NULL
#N3_up$shared.name <- NULL

#renaming column name to GeneID

#colnames(N1_down)[4] = 'GeneID' 
#colnames(N3_down)[4] = 'GeneID' 
#colnames(N1_up)[4] = 'GeneID' 
#colnames(N3_up)[4] = 'GeneID' 

#colnames(N1_down)[5] = 'P-adj' 
#colnames(N3_down)[5] = 'P-adj' 
#colnames(N1_up)[5] = 'P-adj' 
#colnames(N3_up)[5] = 'P-adj' 

#Visualization

#ggmaplot(DEG_NANOS1, main = expression("1" %->% "2"),
#         fdr = 0.01, fc = 0.5, size = 0.4,
#         palette = c("#B31B21", "#1465AC", "darkgray"),
#         genenames = as.vector(DEG_NANOS1$GeneID),
#         legend = "top", top = 20,
#         font.label = c("bold", 11),
#         font.legend = "bold",
#         font.main = "bold",
#         ggtheme = ggplot2::theme_minimal())


#Error in ggmaplot(DEG_NANOS1, main = expression("Group 1" %->% "Group 2"),  : 
#The colnames of data must contain: baseMean, log2FoldChange, padj
#I will change the colnames tomorrow and we will see

#colnames(DEG_NANOS1)[3] = 'log2FoldChange'
#colnames(DEG_NANOS1)[7] = 'padj'
#colnames(DEG_NANOS3)[3] = 'log2FoldChange'
#colnames(DEG_NANOS3)[7] = 'padj'

#Error in ggmaplot(DEG_NANOS1, main = expression("Group 1" %->% "Group 2"),  : 
#genenames should be of length nrow(data). 
#I think I need to clean up NA matrices
#I dont do it anymore the graph looks unproper

#na.omit(DEG_NANOS1) -> DEG_NANOS1
#na.omit(DEG_NANOS3) -> DEG_NANOS3

#omitting NANOS1 and NANOS3 from the dataset so they dont change the scale of the plot
#DEG_NANOS1 <- DEG_NANOS1[-c(1), ]
#DEG_NANOS3 <- DEG_NANOS3[-c(1), ]

#Mapping geneIDs to gene symbols
#So far I can't do this in R version 3.6.3 due to bioMart not supporting it. 
#I managed it. mapping will work now

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")
#library(org.Hs.eg.db)
#?org.Hs.eg

#I miserably failed mapping using R. I might go back to this later
#I will export the datasets as tables import them into cytoscape map there and import back to R

#head(DEG_NANOS1)
#write.table(DEG_NANOS1, file = 'DEG_NANOS1.csv', sep = ',', dec = '.') 

#head(DEG_NANOS3)
#write.table(DEG_NANOS3, file = 'DEG_NANOS3.csv', sep = ',', dec = '.') 

#after this i made a simple manipulation in notepad and added x as the first column because the ranking was
#exported for some reason. i need to solve this sometime

#importing the data back from cytoscape after mapping

DEG_NANOS1_HGNC <- read.table(file = 'DEG_NANOS1_HGNC.csv', header = TRUE, sep = ',', dec = '.')
DEG_NANOS3_HGNC <- read.table(file = 'DEG_NANOS3_HGNC.csv', header = TRUE, sep = ',', dec = '.')
#head(DEG_NANOS1_HGNC) 

#I need to process this once again and remove collumns etc 
DEG_NANOS1_HGNC$SUID <- NULL
DEG_NANOS1_HGNC$selected <- NULL
DEG_NANOS1_HGNC$shared.name <- NULL
DEG_NANOS1_HGNC$x <- NULL
DEG_NANOS3_HGNC$SUID <- NULL
DEG_NANOS3_HGNC$selected <- NULL
DEG_NANOS3_HGNC$shared.name <- NULL
DEG_NANOS3_HGNC$x <- NULL

#GGMAPLOT
#modified for etiuda presentation. export size 450x420
ggmaplot(DEG_NANOS1_HGNC, main = expression("NANOS1 overexpression"),
         fdr = 0.01, fc = 0.5, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(DEG_NANOS1_HGNC$HGNC),
         legend = "top", top = 0,
         font.label = c("bold", 14),
         label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal(base_size = 14))
?ggmaplot
#same as above. modified for etiuda presentation labeled genes removed. export size 460x420
ggmaplot(DEG_NANOS3_HGNC, main = expression("NANOS3 overexpression"),
         fdr = 0.01, fc = 0.5, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(DEG_NANOS3_HGNC$HGNC),
         legend = "top", top = 0,
         font.label = c("bold", 14),
         label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal(base_size = 14))


glimpse(DEG_NANOS3_HGNC)

#to visualize only genes of interests. some reason it didnt work.



#i managed to do it by manually manipulating the adj p values of the genes i want to show and then top12 ing them
DEG_NANOS3_HGNC_kamila <- read_xlsx(path = 'C:/Users/Erkut Ilaslan/Desktop/Book1.xlsx', col_names = TRUE)
DEG_NANOS3_HGNC_kamila$SUID <- NULL
DEG_NANOS3_HGNC_kamila$selected <- NULL
DEG_NANOS3_HGNC_kamila$`shared name` <- NULL
DEG_NANOS3_HGNC_kamila$x <- NULL

ggmaplot(DEG_NANOS3_HGNC_kamila, main = expression("NANOS3 OE RNA-Seq"),
         fdr = 0.05, fc = 0.5, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(DEG_NANOS3_HGNC_kamila$HGNC),
         legend = "top", top = 8,
        #label.select = c("BUB1", "CD83"), 
         font.label = c("bold", 18),
         label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal(base_size = 20)
      )
ggsave()
?ggmaplot

#521x563 are the dimensions
