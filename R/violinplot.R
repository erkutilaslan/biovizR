#' Violin plot to visualize biological data.
#'
#' This function visualize biological data as a violin plot
#'
#' @param violin_data Path to the input file.
#' @param x_lab Default is "conditions". Parameter to set label for conditions.
#' @param y_lab Default is "values". Parameter to set label for values.
#' @return A violin plot of the input data.
#' @export


violinplot <- function(violin_data, x_lab = "conditions", y_lab = "values") {

#data import
if (is.character(violin_data) == TRUE) {

violin_data <- readxl::read_xlsx(violin_data, sheet = 1)

}

#data wrangling.
col_names <- colnames(violin_data)
col_length <- length(colnames(violin_data))

violin_data <- tidyr::pivot_longer(violin_data,
				   cols = col_names[1]:col_names[col_length],
				   names_to = x_lab, values_to = y_lab)

#data visuzalization

final_plot <- ggpubr::ggviolin(violin_data,
                               x = x_lab,
                               y = y_lab,
			       xlab = x_lab,
			       ylab = y_lab,
                               fill = x_lab,
                               palette = "jco",
                               add = "boxplot",
                               add.params = list(fill = "white")) #+
#              ggpubr::stat_compare_means(comparisons = "Cell type",
#	                                 label = "p.signif") + #Add significance levels
#              ggpubr::stat_compare_means(label.y = 50) # add global p-vale

plot(final_plot)
return(final_plot)

}

