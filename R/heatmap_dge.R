#' Heatmap for differential gene expression data.
#'
#' This function visualizes differential gene expresion data
#' as heatmap.
#'
#' @param heatmap_data Path to the input file.
#' @param gene_list Path to the list of genes file which will be visualized.
#' @param header Title for the heatmap.
#' @return A heatmap of heatmap_data.
#' @import ComplexHeatmap
#' @export

heatmap_dge <- function(heatmap_data,
                        gene_list = "",
                        header = "") {

  # data import
  if (is.character(heatmap_data) == TRUE) {
    if (grepl(".csv", as.character(heatmap_data)) == TRUE) {
      heatmap_data <- read.csv(heatmap_data, header = TRUE)
    } else {
      heatmap_data <- readxl::read_excel(heatmap_data, sheet = 1)
    }
  }

  # NA omit otherwise we can't move HGNC to rownames
  heatmap_data <- na.omit(heatmap_data)

  # removing row names
  heatmap_data <- tibble::remove_rownames(heatmap_data)

  # moving HGNC to rownames
  colnames(heatmap_data)[1] <- "hgnc_symbol"
  heatmap_data <- tibble::column_to_rownames(heatmap_data, var = "hgnc_symbol")

  # z-score calculation
  heatmap_data <- t(scale(t(heatmap_data)))

  # selecing of genes to visualize e.g. infertility vectors
  if (gene_list != "") {
    gene_list <- read.table(gene_list, header = FALSE, sep = " ")
    gene_list <- gene_list$V1


    # need to take the columns with the genes i want to visualize
    heatmap_data <- heatmap_data[as.character(gene_list), ]
  }

  # removing na values once a mfing time
  heatmap_data <- na.omit(heatmap_data)

  # Visualization
  ComplexHeatmap::Heatmap(heatmap_data,
    heatmap_legend_param = list(
      title = "z-score",
      legend_height = grid::unit(5, "cm"),
      grid_width = grid::unit(0.75, "cm"),
      labels_gp = grid::gpar(fontsize = 12),
      title_gp = grid::gpar(fontsize = 14)
    ),
    column_title_gp = grid::gpar(fontsize = 24),
    column_title = header,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 14),
    column_names_gp = grid::gpar(fontsize = 16),
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    clustering_method_rows = "ward.D2",
    clustering_distance_rows = "pearson"
  )
}
