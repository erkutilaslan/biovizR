#' Analyze and visualize qPCR results with barplot.
#'
#' This function visualizes qPCR results as bar plot.
#'
#' @param qpcr_data Input qPCR data to visualize.
#' @param top Turned off by default. Parameter to set top amount of processes to be visualized.
#' @param go_process Turned off by default. Parameter to add a specific keyword to only visuzalize GO terms that contains it.
#' @return A bar plot of gene ontology results.
#' @import tidyverse
#' @export

barplot_qpcr <- function(qpcr_data,
			 type = "biorad",
			 ref = "",
			 target = "",
			 stat = "t-test") {

#data import 
qpcr_data <- readxl::read_xlsx(qpcr_data, sheet = 1)

#data process

#statistical analysis

#visualization
final_plot <- ggpubr::ggbarplot(qpcr_data,
          x = "GOTerm",
          y = "log10_padj",
          fill = "darkgray",
          xlab = "GO Term",
	  ylab = "p-adjusted(-log10)",
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          label = go_data$`Nr. Genes`,
          title = "",
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          rotate = FALSE,
          ggtheme = ggpubr::theme_pubr(base_size = 18))

plot(final_plot)
return(final_plot)

}
