#' Violin plot to visualize biological data.
#'
#' This function visualize biological data as a violin plot
#'
#' @param violin_data Path to the input file.
#' @return A violin plot of the input data.
#' @export

library(tidyverse)
library(ggpubr)

violinplot <- function(violin_data) {

#data import
violin_data <- readxl::read_xlsx(violin_data, sheet = 1)
violin_data <- readxl::read_xlsx("~/downloads/Toshi_3'UTR.xlsx", sheet = 1)

#data wrangling. try positional arrangement will be better for function

violin_data <- tidyr::pivot_longer(violin_data, cols = hPSC.SB:hPGC_Day4, names_to = "Cell type", values_to = "UTR_lenght")

#data visuzalization

final_plot <- ggpubr::ggviolin(test,
                               x = "Cell type",
                               y = "UTR_lenght",
			       ylab = "UTR lenght",
                               fill = "Cell type",
                               palette = "jco",
                               add = "boxplot",
                               add.params = list(fill = "white")) #+
#              ggpubr::stat_compare_means(comparisons = "Cell type",
#	                                 label = "p.signif") + #Add significance levels
#              ggpubr::stat_compare_means(label.y = 50) # add global p-vale

plot(final_plot)
return(final_plot)


}

