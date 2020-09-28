devtools::install_github("erkutilaslan/biovizR")

library(biovizR)

barplot_go("C://Users/Anelor/Desktop/go_data.xls", top = 10, go_process = "cell cycle")
#works as intended

maplot_dge("C://Users/Anelor/Desktop/maplot_data.csv", FDR = 0.05, FC= 1)
#works as intended maybe I can add gene names here.

heatmap_dge(heatmap_data = "C://Users/Anelor/Desktop/heatmap_data.csv", List = "C://Users/Anelor/Desktop/test_File_heatmap.txt")
#works as intended its just perfect.

bubbleplot_go("C://Users/Anelor/Desktop/go_data.xls", top = 10)
#Error in scale_size(range = c(2, 12)) : 

dotplot_ipms(ip_data = "C://Users/Anelor/Desktop/ipms_data.xlsx")
#Error in aes(label = ifelse(avg.ctrl < 2, as.character(HGNC), "")) :

remove.packages("biovizR")
