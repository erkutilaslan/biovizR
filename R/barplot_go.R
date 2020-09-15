#' Bar plot for gene ontology visualization.
#'
#' This function visualizes gene ontology results from ClueGO as bar plot.
#'
#' @param go_data Input GO data to visualize.
#' @param top Turned off by default. Parameter to set top amount of processes to be visualized.
#' @param go_process Turned off by default. Parameter to add a specific keyword to only visuzalize GO terms that contains it.
#' @return A bar plot of gene ontology results.
#' @import tidyverse
#' @export


barplot_go <- function(go_data, top = " ", go_process = " ") {

# import data
go_data <- readxl::read_excel(go_data, sheet = 1)

# selecting GOTerms for visualization
if (go_process != " ") {
go_data <- go_data[grep(go_process,
                       go_data$GOTerm,
                       ignore.case = TRUE), ]
}

# deduplication. because of some terms are present in multiple
# GO groups they are duplicated

# first removing 'group' columns they interfere with deduplication

go_data <- dplyr::select(go_data, -6, -7, -8, -9)

# deduplication

go_data <- dplyr::distinct(go_data)

# before visualization I need to convert adj-p by -log10 for visualization

# first converting column name so its easier to mutate

colnames(go_data)[5] <- "padj"

# -log10 conversion of padj by mutate

go_data <- dplyr::mutate(log10_padj = -log10(padj))

# for long GO lists I want to visualize only top 10/20.
# this is not possible with top = 10 command
# while generating plots. so i will generate new tables of top 10 top 20.

go_data <- dplyr::arrange(go_data, dplyr::desc(go_data$log10_padj))

if (top != " ") {

go_data <- dplyr::slice(go_data, 1:top)

}

 #visualization
final_plot <- ggpubr::ggbarplot(go_data,
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
          rotate = TRUE,
          ggtheme = ggpubr::theme_pubr(base_size = 18))

plot(final_plot)
return(final_plot)

}
