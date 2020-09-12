#' MA-plot visualization of RNA-Seq DGE analysis.
#'
#' This function visualizes RNA-Seq DGE results as an MA-plot.
#'
#' @param rnaseq_data Path to the input file.
#' @param FDR Default 0.01. Adjust FDR value.
#' @param FC Default 0.5. Adjust log2FC threshold.
#' @return An ma-plot of RNA-Seq DGE data.
#' @export

maplot_rnaseq <- function(rnaseq_data, FDR = 0.01, FC = 0.5) {

#data import
rnaseq_data <- read.table(file = rnaseq_data,
			  header = TRUE, sep = ",", dec = ".")

#head(dat)

#I need to process this once again and remove collumns etc 
rnaseq_data$SUID <- NULL
rnaseq_data$selected <- NULL
rnaseq_data$shared.name <- NULL
rnaseq_data$x <- NULL

# Visuzalization

final_plot <- ggpubr::ggmaplot(rnaseq_data, main = " ",
         fdr = FDR, fc = FC, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(rnaseq_data$HGNC),
         legend = "top", top = 0,
         font.label = c("bold", 14),
         label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal(base_size = 14))

plot(final_plot)
return(final_plot)

}
