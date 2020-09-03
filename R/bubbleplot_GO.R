library(tidyverse)
library(readxl)
library(ggthemes)
library(magrittr)
library(ggpubr)
library(biomaRt)
library(clusterProfiler)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#	install.packages("BiocManager")
#BiocManager::install("biomaRt")

#GO looking for the same terms in every clone
#also visualizing the ones she liked the most
#weekend project

#well i decided to find the genes differentially expressed

setwd('C:/Users/Erkut Ilaslan/Desktop/Data Integration/Karolina Rassek/RNAseq/')

#importing data

#JpC <- read.delim('JpC.txt', dec = ',')
JpC11 <- read.delim('JpC11.txt', dec = ',')
JpC13 <- read.delim('JpC13.txt', dec = ',')
JpC14 <- read.delim('JpC14.txt', dec = ',')
JpC23 <- read.delim('JpC23.txt', dec = ',')
JpC26 <- read.delim('JpC26.txt', dec = ',')
JpP13 <- read.delim('JpP13.txt', dec = ',')
JpP15 <- read.delim('JpP15.txt', dec = ',')
JpP22 <- read.delim('JpP22.txt', dec = ',')
JpP23 <- read.delim('JpP23.txt', dec = ',')
JpP25 <- read.delim('JpP25.txt', dec = ',')

#need to remove Description column because it is NA in all them and
#it interferes with na.omit after left_join

JpC11 <- JpC11[ , -9]
JpC13 <- JpC13[ , -9]
JpC14 <- JpC14[ , -9]
JpC23 <- JpC23[ , -9]
JpC26 <- JpC26[ , -9]

JpP13 <- JpP13[ , -9]
JpP22 <- JpP22[ , -9]
JpP23 <- JpP23[ , -9]
JpP25 <- JpP25[ , -9]


#left join
JpC <- left_join(JpC11, JpC13, 
                 by = c('Gene' = 'Gene'),
                 copy = FALSE)

JpC <- left_join(JpC, JpC14, 
                 by = c('Gene' = 'Gene'),
                 copy = FALSE)

JpC <- left_join(JpC, JpC23, 
                 by = c('Gene' = 'Gene'),
                 copy = FALSE)

#JpC <- left_join(JpC, JpC26, 
#                 by = c('Gene' = 'Gene'),
#                 copy = FALSE)


####only best 2 
JpC <- left_join(JpC13, JpC23, 
                 by = c('Gene' = 'Gene'),
                 copy = FALSE)

#################

JpP <- left_join(JpP13, JpP22,
                 by = c('Gene' = 'Gene'),
                 copy = FALSE)

#JpP <- left_join(JpP, JpP23,
#                 by = c('Gene' = 'Gene'),
#                 copy = FALSE)

JpP <- left_join(JpP, JpP25,
                 by = c('Gene' = 'Gene'),
                 copy = FALSE)

#na.omit

JpC <- na.omit(JpC)

JpP <- na.omit(JpP)

#export for GO in Cytoscape using ClueGO

write.csv(JpC, file = 'JpC_genes_only3best.csv', quote = FALSE)

write.csv(JpP, file = 'JpP_genes_only3best.csv', quote = FALSE)

#importing GO results for visualization

JpC_GO <- read_xls('JpC_GO_only2best.xls')

#i need to copy the script from mine and continue with visualization.

#https://rpubs.com/Koundy/71792
#https://github.com/jrnold/ggthemes

#processing data for visualization

JpC_GO %>% 
  dplyr::select(-6,-7,-8,-9) %>% #removing group named columns
  distinct() -> x #deduplication
colnames(x)[5] = 'padj'
colnames(x)[6] = 'associated_genes_percent'
colnames(x)[7] = 'nr.genes' #changing colname for easier mutate
x %>%
  mutate('log10_padj' = -log10(padj)) -> x #calculating -log10 value 

#i wanted to try how to  sort for bubble plots but it didnt work
#x <- arrange(x, x$padj)

ggbarplot(x, 
          x = 'GOTerm', 
          y = 'log10_padj', 
          fill = "#00BFFF",            
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          label = 'nr.genes',
          title = 'TMEM 244 Overexpression',
          lab.size = 5,
          lab.vjust = 0.4,
          lab.hjust = 1.2,
         # title(main = ),
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          rotate = TRUE,
          ggtheme = theme_pubr(base_size = 18))

??ggbarplot

#I will try bubbleplot
#https://stackoverflow.com/questions/26757026/bubble-chart-with-ggplot2
#https://www.r-graph-gallery.com/320-the-basis-of-bubble-plot.html
#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/


ggplot(x, aes(x = log10_padj, y = GOTerm, label = nr.genes)) +
  geom_point(aes(size = nr.genes)) + 
  #geom_text(hjust = 1, size = 2) +
  scale_size(range = c(1,6)) +
  theme_bw()
  
#lets try ggpubr

ggscatter(x, x = 'log10_padj', y = 'GOTerm',
          size = "nr.genes",
          color = '#00AFBB'
          ) +
  scale_size(range = c(2, 12))    # Adjust the range of points size
  
?ggscatter
