#' Heatmap for differential gene expression data.
#'
#' This function visualizes differential gene expresion data
#' as heatmap.
#'
#' @param heatmap_data Path to the input file.
#' @param List Path to the list of genes file which will be visualized.
#' @return A heatmap of heatmap_data.
#' @import ComplexHeatmap
#' @export


heatmap_dge <- function(heatmap_data, List = "") {

#data import

heatmap_data <- read.csv(file = heatmap_data, header = TRUE)

#NA omit otherwise we can't move HGNC to rownames

heatmap_data  <- na.omit(heatmap_data)

#removing row names

heatmap_data <- tibble::remove_rownames(heatmap_data)

#moving HGNC to rownames
#this is an issue point here. because we need to able to move HGNC to the first row.
#probably i should name the first row into something automatically as HGNC and it should
#accept such data for now and later on we will make sure that it will accept
# all kind of identifiersssss!
heatmap_data <- tibble::column_to_rownames(heatmap_data, var = "hgnc_symbol")

#z-score calculation

heatmap_data <- t(scale(t(heatmap_data)))

#selecing of genes to visualize e.g. infertility vectors
if (List != "") {

List <- read.table(List, header = FALSE, sep = " ")

List <- List$V1

#need to take the columns with the genes i want to visualize

heatmap_data <- heatmap_data[as.character(List), ]

}

#removing na values once a mfing time
heatmap_data  <- na.omit(heatmap_data)
# Visualization

final_plot <- ComplexHeatmap::Heatmap(
  heatmap_data,
  #column_title = "NANOS1: Male Infertility",
  heatmap_legend_param = list(title = "z-score",
                              legend_height = grid::unit(5, "cm"),
                              grid_width = grid::unit(0.75, "cm"),
                              labels_gp = grid::gpar(fontsize = 12),
                              title_gp = grid::gpar(fontsize = 14)),
  column_title_gp = grid::gpar(fontsize = 18),
  #col = greenred(75),
  #show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 14),
  column_names_gp = grid::gpar(fontsize = 16),
  #show_column_names = TRUE
  #cluster_rows = FALSE,
  cluster_columns = FALSE,
  #show_column_dend = TRUE,
  show_row_dend = FALSE,
  #row_dend_reorder = TRUE,
  #column_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = "pearson"
  #clustering_method_columns = "ward.D2",
  #width = unit(150, "mm"),
  #height = unit(130, "mm")
)

return(final_plot)

}
