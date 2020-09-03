library(ggthemes)
library(magrittr)
library(ggpubr)
library(biomaRt)
library(readxl)
library(tidyverse)

#GGBARPLOTS for GO results

#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/

#first import data

N1_down_GO <- read_excel('~/nanos_test.xlsx', sheet = 1)
N1_up_GO <- read_excel('N1_up_GO.xls', sheet = 1)
N3_down_GO <- read_excel('N3_down_GO.xls', sheet = 1)
N3_up_GO <- read_excel('N3_up_GO.xls', sheet = 1)

#selecting GOTerms for visualization

N1_down_GO_cellcycle <- N1_down_GO[grep('cell cycle', 
                                        N1_down_GO$GOTerm, 
                                        ignore.case = TRUE), ]

N1_down_GO_development <- N1_down_GO[grep('development', 
                                          N1_down_GO$GOTerm, 
                                          ignore.case = TRUE), ]


N1_down_GO_differentiation <- N1_down_GO[grep('differentiation', 
                                          N1_down_GO$GOTerm, 
                                          ignore.case = TRUE), ]

N3_down_GO_cellcycle <- N3_down_GO[grep('cell cycle', 
                                        N3_down_GO$GOTerm, 
                                        ignore.case = TRUE), ]

N3_down_GO_development <- N3_down_GO[grep('development', 
                                          N3_down_GO$GOTerm, 
                                          ignore.case = TRUE), ]

N3_down_GO_differentiation <- N3_down_GO[grep('differentiation', 
                                              N3_down_GO$GOTerm, 
                                              ignore.case = TRUE), ]

N3_up_GO_cellcycle <- N3_up_GO[grep('cell cycle', 
                                    N3_up_GO$GOTerm, 
                                    ignore.case = TRUE), ]

N3_up_GO_development <- N3_up_GO[grep('development', 
                                      N3_up_GO$GOTerm, 
                                      ignore.case = TRUE), ]

N3_up_GO_differentiation <- N3_up_GO[grep('differentiation', 
                                      N3_up_GO$GOTerm, 
                                      ignore.case = TRUE), ]

#removing N3_up_GO_cellcycle because there are no GOterms

remove(N3_up_GO_cellcycle)

#optional script for multiple keyword selection using grepl
#N1_GO_test_cellcycle_proliferation <- N1_GO_test[grepl('cell cycle|proliferation', 
#                                                      N1_GO_test$GOTerm, 
#                                                      ignore.case = TRUE), ]



#deduplication. because of some terms are present in multiple GO groups they are duplicated

#first removing 'group' columns they interfere with deduplication

colnames(N1_down_GO_cellcycle) 

N1_down_GO_cellcycle_dedup <- dplyr::select(N1_down_GO_cellcycle, -6,-7,-8,-9)

N1_down_GO_development_dedup <- dplyr::select(N1_down_GO_development, -6,-7,-8,-9)

N1_down_GO_differentiation_dedup <- dplyr::select(N1_down_GO_differentiation, -6,-7,-8,-9)

N3_down_GO_cellcycle_dedup <- dplyr::select(N3_down_GO_cellcycle, -6,-7,-8,-9)

N3_down_GO_development_dedup <- dplyr::select(N3_down_GO_development, -6,-7,-8,-9)

N3_down_GO_differentiation_dedup <- dplyr::select(N3_down_GO_differentiation, -6,-7,-8,-9)

N3_up_GO_development_dedup <- dplyr::select(N3_up_GO_development, -6,-7,-8,-9)

N3_up_GO_differentiation_dedup <- dplyr::select(N3_up_GO_differentiation, -6,-7,-8,-9)

colnames(N1_GO_test_cellcycle_dedup) 

#deduplication
          
N1_down_GO_cellcycle_dedup <- distinct(N1_down_GO_cellcycle_dedup)

N1_down_GO_development_dedup <- distinct(N1_down_GO_development_dedup)

N1_down_GO_differentiation_dedup <- distinct(N1_down_GO_differentiation_dedup)

N3_down_GO_cellcycle_dedup <- distinct(N3_down_GO_cellcycle_dedup)

N3_down_GO_development_dedup <- distinct(N3_down_GO_development_dedup)

N3_down_GO_differentiation_dedup <- distinct(N3_down_GO_differentiation_dedup)

N3_up_GO_development_dedup <- distinct(N3_up_GO_development_dedup)

N3_up_GO_differentiation_dedup <- distinct(N3_up_GO_differentiation_dedup)

#before visualization I need to convert adj-p by -log10 for visualization

#first converting column name so its easier to mutate

colnames(N1_down_GO_cellcycle_dedup)[5] = 'padj'

colnames(N1_down_GO_development_dedup)[5] = 'padj'

colnames(N1_down_GO_differentiation_dedup)[5] = 'padj'

colnames(N3_down_GO_cellcycle_dedup)[5] = 'padj'

colnames(N3_down_GO_development_dedup)[5] = 'padj'

colnames(N3_down_GO_differentiation_dedup)[5] = 'padj'

colnames(N3_up_GO_development_dedup)[5] = 'padj'

colnames(N3_up_GO_differentiation_dedup)[5] = 'padj'

# -log10 conversion of padj by mutate

N1_down_GO_cellcycle_dedup %>% 
  mutate('log10_padj' = -log10(padj)) -> N1_down_GO_cellcycle_dedup

N1_down_GO_development_dedup %>% 
  mutate('log10_padj' = -log10(padj)) -> N1_down_GO_development_dedup

N1_down_GO_differentiation_dedup %>% 
  mutate('log10_padj' = -log10(padj)) -> N1_down_GO_differentiation_dedup

N3_down_GO_cellcycle_dedup %>% 
  mutate('log10_padj' = -log10(padj)) -> N3_down_GO_cellcycle_dedup

N3_down_GO_development_dedup %>% 
  mutate('log10_padj' = -log10(padj)) -> N3_down_GO_development_dedup

N3_down_GO_differentiation_dedup %>% 
  mutate('log10_padj' = -log10(padj)) -> N3_down_GO_differentiation_dedup

N3_up_GO_development_dedup %>% 
  mutate('log10_padj' = -log10(padj)) -> N3_up_GO_development_dedup

N3_up_GO_differentiation_dedup %>% 
  mutate('log10_padj' = -log10(padj)) -> N3_up_GO_differentiation_dedup

#for long GO lists I want to visualize only top 10/20. this is not possible with top = 10 command
#while generating plots. so i will generate new tables of top 10 top 20.

N1_down_GO_development_dedup %>% 
  arrange(desc(log10_padj)) %>%
  slice(1:10) -> N1_down_GO_development_dedup_top10

N3_down_GO_development_dedup %>% 
  arrange(desc(log10_padj)) %>%
  slice(1:10) -> N3_down_GO_development_dedup_top10

N3_up_GO_development_dedup %>% 
  arrange(desc(log10_padj)) %>%
  slice(1:10) -> N3_up_GO_development_dedup_top10

#adding enter to long terms by \n in between

N1_down_GO_cellcycle_dedup[11,2] = 'DNA damage response, signal transduction by p53 \n class mediator resulting in cell cycle arrest'

N1_down_GO_cellcycle_dedup[9,2] = 'signal transduction involved in mitotic cell cycle \n checkpoint'


#visualization

?ggbarplot()
?theme_pubr



ggbarplot(N1_down_GO_cellcycle_dedup, 
          x = 'GOTerm', 
          y = 'log10_padj', 
          fill = "darkgray",            
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          label = N1_down_GO_cellcycle_dedup$`Nr. Genes`,
          title = 'NANOS1 Downregulated GO:Cell cycle',
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          rotate = TRUE,
          ggtheme = theme_pubr(base_size = 18))

ggbarplot(N1_down_GO_development_dedup_top10, x = 'GOTerm', y = 'log10_padj', 
          fill = "darkgray",            
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          label = N1_down_GO_development_dedup_top10$`Nr. Genes`,
          title = 'Top 10 NANOS1 Downregulated GO:Development',
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          position.title  = -10,
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          rotate = TRUE,
          ggtheme = theme_pubr(base_size = 18)) #size of text

ggbarplot(N1_down_GO_differentiation_dedup, x = 'GOTerm', y = 'log10_padj', 
          fill = "darkgray",            
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          label = N1_down_GO_differentiation_dedup$`Nr. Genes`,
          title = 'NANOS1 Downregulated GO:Differentiation',
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          position.title  = -10,
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          rotate = TRUE,
          ggtheme = theme_pubr(base_size = 18)) #size of text

ggbarplot(N3_down_GO_cellcycle_dedup, x = 'GOTerm', y = 'log10_padj', 
          fill = "darkgray",            
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          label = N3_down_GO_cellcycle_dedup$`Nr. Genes`,
          title = 'NANOS3 Downregulated GO:Cell cycle',
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          rotate = TRUE,
          ggtheme = theme_pubr(base_size = 18) )

ggbarplot(N3_down_GO_differentiation_dedup, x = 'GOTerm', y = 'log10_padj', 
          fill = "darkgray",            
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          label = N3_down_GO_differentiation_dedup$`Nr. Genes`,
          title = 'NANOS3 Downregulated GO:Differentiation',
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          rotate = TRUE,
          ggtheme = theme_pubr(base_size = 18) )

ggbarplot(N3_up_GO_development_dedup_top10, x = 'GOTerm', y = 'log10_padj', 
          fill = "darkgray",            
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          label = N3_up_GO_development_dedup_top10$`Nr. Genes`,
          title = 'Top 10 NANOS3 Upregulated GO:Development',
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          rotate = TRUE,
          ggtheme = theme_pubr(base_size = 18) )

ggbarplot(N3_up_GO_differentiation_dedup, x = 'GOTerm', y = 'log10_padj', 
          fill = "darkgray",            
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          label = N3_up_GO_differentiation_dedup$`Nr. Genes`,
          title = 'NANOS3 Upregulated GO:Differentiation',
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          rotate = TRUE,
          ggtheme = theme_pubr(base_size = 18) )

# jan et al., 2017 process this and add to germ cell specific genes database and check for clusters once again

