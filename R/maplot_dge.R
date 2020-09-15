#' MA-plot visualization of DGE analysis.
#'
#' This function visualizes DGE results as an MA-plot.
#'
#' @param dge_data Path to the input file.
#' @param FDR Default 0.01. Adjust FDR value.
#' @param FC Default 0.5. Adjust log2FC threshold.
#' @return An ma-plot of RNA-Seq DGE data.
#' @export

maplot_dge <- function(dge_data, FDR = 0.01, FC = 0.5) {

#data import
dge_data <- read.table(file = dge_data,
			  header = TRUE, sep = ",", dec = ".")

#head(dat)

#I need to process this once again and remove collumns etc 
dge_data$SUID <- NULL
dge_data$selected <- NULL
dge_data$shared.name <- NULL
dge_data$x <- NULL

# Visuzalization

final_plot <- ggpubr::ggmaplot(dge_data, main = " ",
         fdr = FDR, fc = (FC), size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(dge_data$HGNC),
         legend = "top", top = 0,
         font.label = c("bold", 14),
         label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal(base_size = 14))

plot(final_plot)
return(final_plot)

}
