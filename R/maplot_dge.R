#' MA-plot visualization of DGE analysis.
#'
#' This function visualizes DGE results as an MA-plot.
#'
#' @param dge_data Path to the input file.
#' @param FDR Default 0.01. Adjust FDR value.
#' @param FC Default 0.5. Adjust log2FC threshold.
#' @param TOP Default 10. Adjust top number of DE genes to visualize.
#' @param header Default is empty. Set a title for the MA-plot.
#' @param type Default deseq2. Specify input datatype.
#' @return An ma-plot of RNA-Seq DGE data.
#' @export


maplot_dge <- function(dge_data, FDR = 0.01, FC = 0.5, TOP = 10, type = "deseq2", header = "") {

#data import
if (is.character(dge_data) == TRUE) {
 
  if (grepl(".csv", as.character(dge_data)) == TRUE) {
    
    dge_data <- read.csv(dge_data) 

  } else {

    dge_data <- readxl::read_excel(dge_data, sheet = 1)

  }

}

if (type == "deseq2") {

  #processing colnames for visualization
  colnames(dge_data)[2] = "baseMean"
  colnames(dge_data)[3] = "log2FoldChange"
  colnames(dge_data)[7] = "padj"

}

if (type == "edger") {



}

# Visuzalization
final_plot <- ggpubr::ggmaplot(dge_data, main = header,
         fdr = FDR, fc = FC, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(dge_data$HGNC),
         legend = "top",
	 top = TOP,
         font.label = c("bold", 14),
         label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal(base_size = 14))

plot(final_plot)
return(final_plot)

}
