#' Bar plot for gene ontology visualization.
#'
#' This function visualizes gene ontology results from ClueGO as bar plot.
#'
#' @param go_data Input GO data to visualize.
#' @param type Default cluego. Parameter to specify input data source.
#' @param top Default value is 50. Parameter to set top number of processes to be visualized.
#' @param go_process Turned off by default. Parameter to add a specific keyword to only visuzalize GO terms that contains it.
#' @param min_genes Turned off by default. Parameter to set threshold for biological processes containing minimum number of genes.
#' @param header Set a title for the barplot.
#' @return A bar plot of gene ontology results.
#' @import tidyverse
#' @export

barplot_go <- function(go_data,
		       type = "cluego",
		       top = "50",
		       go_process = "",
		       min_genes = "",
		       header = "") {

# import data
if (is.character(go_data) == TRUE) {
 
	if (grepl(".csv", as.character(go_data)) == TRUE) {
    
		go_data <- read.csv(go_data) 

	} else {

		go_data <- readxl::read_excel(go_data, sheet = 1)

	}

}

if (type == "cluego") {

	#removing unnecessary columns they interfere with deduplication
	go_data <- dplyr::select(go_data, -6, -7, -8, -9)

	# deduplication
	go_data <- dplyr::distinct(go_data)

	# converting column name so its easier to mutate
	colnames(go_data)[5] <- "padj"

	# -log10 conversion of padj by mutate
	go_data <- dplyr::mutate(go_data, log10_padj = -log10(go_data$padj))

	# arranging the GO Terms in ascending order
	go_data <- dplyr::arrange(go_data, log10_padj)

}

if (type == "panther") {

	#removing unnecessary columns they interfere with deduplication
	go_data <- dplyr::select(go_data, 1, 3, 8)

	# deduplication
	go_data <- dplyr::distinct(go_data)

	# converting column name so its easier to mutate
	colnames(go_data)[c(1, 2, 3)] <- c("GOTerm", "Nr. Genes", "padj")

	# -log10 conversion of padj by mutate
	go_data <- dplyr::mutate(go_data, log10_padj = -log10(go_data$padj))

	# arranging the GO Terms in ascending order
	go_data <- dplyr::arrange(go_data)

}

if (type == "david") {

	# first removing unnecessary columns they interfere with deduplication
	go_data <- dplyr::select(go_data, 1, 2, 3, 13)
  
	# deduplication
	go_data <- dplyr::distinct(go_data)

	# converting column names
	colnames(go_data)[c(2, 3, 4)] <- c("GOTerm", "Nr. Genes", "padj")

	# -log10 conversion of padj by mutate
	go_data <- dplyr::mutate(go_data, log10_padj = -log10(go_data$padj))

	# arranging the GO Terms in ascending order
	go_data <- dplyr::arrange(go_data)

}

if (type == "gprofiler") {

	#removing unnecessary columns they interfere with deduplication
	go_data <- dplyr::select(go_data, 3, 2, 1, 5, 8)

	#renaming columns for visualization
	colnames(go_data)[c(2,4,5)] <- c("GOTerm", "log10_padj", "Nr. Genes")

	# deduplication
	go_data <- dplyr::distinct(go_data)

}

# selecting GOTerms for visualization
if (go_process != " ") {

	go_data <- go_data[grep(go_process,
    	                        go_data$GOTerm,
                                ignore.case = TRUE), ]
}

#selecting top processes for visualization

if (top != " ") {

	go_data <- dplyr::slice(go_data, 1:top)

}

if (min_genes != " ") {

	go_data <- dplyr::filter(go_data, go_data$`Nr. Genes` >= min_genes)

}

#visualization
final_plot <- ggpubr::ggbarplot(go_data,
                                x = "GOTerm",
                                y = "log10_padj",
                                fill = "darkgray",
                                xlab = "GO Term",
                                ylab = "p-adjusted(-log10)",
                                size = 0.5,
                                palette = "jco",
                                label = go_data$`Nr. Genes`,
                                title = header,
                                lab.size = 5,
                                lab.vjust = 0.5,
                                lab.hjust = 1.2,
                                sort.by.groups = FALSE,
                                rotate = TRUE,
                                ggtheme = ggpubr::theme_pubr(base_size = 18))

plot(final_plot)
return(final_plot)

}
