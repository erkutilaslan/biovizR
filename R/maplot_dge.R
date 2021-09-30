#' MA-plot visualization of DGE analysis.
#'
#' This function visualizes DGE results as an MA-plot.
#'
#' @param dge_data Path to the input file.
#' @param FDR Default 0.05. Adjust FDR value.
#' @param FC Default 2. Adjust log2FC threshold.
#' @param TOP Default 10. Adjust top number of DE genes to visualize.
#' @param header Default is empty. Set a title for the MA-plot.
#' @param type Default deseq2. Specify input datatype.
#' @param mapping Default is FALSE. Set to TRUE to convert the gene identifiers to HGNC.
#' @return An ma-plot of RNA-Seq DGE data.
#' @export

dge_data <- "~/transcriptomics/karolina_rassek/rna-seq_new/dge/4vs4/TMEM244_expression_High_vs_Control.tabular"
FDR <- 0.05
FC <- 2
TOP <- 10
mapping <- FALSE
type <- "deseq2"
header <- "biovizrrr"

maplot_dge <- function(dge_data,
		       FDR = 0.05,
		       FC = 2,
		       TOP = 10,
		       type = "deseq2",
		       header = "",
		       mapping = FALSE) {

#mapping function
map_to_hgnc <- function (map_data) {

	if (grepl("^\\d+$", map_data[1, 1]) == TRUE) {

		ensembl <- biomaRt::useMart("ensembl")
		ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

		datasets <- biomaRt::listDatasets(ensembl)
		attributes <- biomaRt::listAttributes(ensembl)
		filters <- biomaRt::listFilters(ensembl)

		# prepare query identifiers from the transcripts dataset
		gene_list <- dge_data[ ,1]

		ID_to_HGNC <- biomaRt::getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
       				            filters = c(filters = "entrezgene_id"),
       			         	    values = gene_list,
					    uniqueRows = TRUE,
			                    useCache = FALSE,
			                    mart = ensembl)

		#removing duplicated rows
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]
	
		# replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
		ID_to_HGNC[ID_to_HGNC == ""] <- NA
		ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC),]

		#add HGNC to dataframe by joining them
		dge_data <- dplyr::left_join(dge_data, ID_to_HGNC, 
       			                     by = c("GeneID" = 'entrezgene_id'),
                            		     copy = FALSE)

		#replace gene_id column with new identifier
		dge_data <- dge_data[, -1]
		dge_data <- dge_data[, c(7,1:6)]

	} else if (grepl("^(\\w+\\d+(\\.\\d+)?)|(NP_\\d+)$", map_data[1,1]) == TRUE) {

		#NCBI Protein
		ensembl <- biomaRt::useMart("ensembl")
		ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

		datasets <- biomaRt::listDatasets(ensembl)
		attributes <- biomaRt::listAttributes(ensembl)
		filters <- biomaRt::listFilters(ensembl)

		# prepare query identifiers from the transcripts dataset
		gene_list <- dge_data[ ,1]

		ID_to_HGNC <- biomaRt::getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
       				            filters = c(filters = "entrezgene_id"),
       			         	    values = gene_list,
					    uniqueRows = TRUE,
			                    useCache = FALSE,
			                    mart = ensembl)

		#removing duplicated rows
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]
	
		# replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
		ID_to_HGNC[ID_to_HGNC == ""] <- NA
		ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC),]

		#add HGNC to dataframe by joining them
		dge_data <- dplyr::left_join(dge_data, ID_to_HGNC, 
       			                     by = c("GeneID" = 'entrezgene_id'),
                            		     copy = FALSE)

		#replace gene_id column with new identifier
		dge_data <- dge_data[, -1]
		dge_data <- dge_data[, c(7,1:6)]



	} else if (grepl("^RF\\d{5}$", map_data[1,1]) == TRUE) {

		#RFAM
		ensembl <- biomaRt::useMart("ensembl")
		ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

		datasets <- biomaRt::listDatasets(ensembl)
		attributes <- biomaRt::listAttributes(ensembl)
		filters <- biomaRt::listFilters(ensembl)

		# prepare query identifiers from the transcripts dataset
		gene_list <- dge_data[ ,1]

		ID_to_HGNC <- biomaRt::getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
       				            filters = c(filters = "entrezgene_id"),
       			         	    values = gene_list,
					    uniqueRows = TRUE,
			                    useCache = FALSE,
			                    mart = ensembl)

		#removing duplicated rows
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]
	
		# replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
		ID_to_HGNC[ID_to_HGNC == ""] <- NA
		ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC),]

		#add HGNC to dataframe by joining them
		dge_data <- dplyr::left_join(dge_data, ID_to_HGNC, 
       			                     by = c("GeneID" = 'entrezgene_id'),
                            		     copy = FALSE)

		#replace gene_id column with new identifier
		dge_data <- dge_data[, -1]
		dge_data <- dge_data[, c(7,1:6)]




	} else if (grepl("^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\\.\\d+)?$", map_data == TRUE)) {

		#UNIPROT
		ensembl <- biomaRt::useMart("ensembl")
		ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

		datasets <- biomaRt::listDatasets(ensembl)
		attributes <- biomaRt::listAttributes(ensembl)
		filters <- biomaRt::listFilters(ensembl)

		# prepare query identifiers from the transcripts dataset
		gene_list <- dge_data[ ,1]

		ID_to_HGNC <- biomaRt::getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
       				            filters = c(filters = "entrezgene_id"),
       			         	    values = gene_list,
					    uniqueRows = TRUE,
			                    useCache = FALSE,
			                    mart = ensembl)

		#removing duplicated rows
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]
	
		# replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
		ID_to_HGNC[ID_to_HGNC == ""] <- NA
		ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC),]

		#add HGNC to dataframe by joining them
		dge_data <- dplyr::left_join(dge_data, ID_to_HGNC, 
       			                     by = c("GeneID" = 'entrezgene_id'),
                            		     copy = FALSE)

		#replace gene_id column with new identifier
		dge_data <- dge_data[, -1]
		dge_data <- dge_data[, c(7,1:6)]




	} else if (grepl("^((ENS[FPTG]\\d{11}(\\.\\d+)?)|(FB\\w{2}\\d{7})|(Y[A-Z]{2}\\d{3}[a-zA-Z](\\-[A-Z])?)|([A-Z_a-z0-9]+(\\.)?(t)?(\\d+)?([a-z])?))$", map_data == TRUE)) {
		
		#ENSEMBL
		ensembl <- biomaRt::useMart("ensembl")
		ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

		datasets <- biomaRt::listDatasets(ensembl)
		attributes <- biomaRt::listAttributes(ensembl)
		filters <- biomaRt::listFilters(ensembl)

		# prepare query identifiers from the transcripts dataset
		gene_list <- dge_data[ ,1]

		ID_to_HGNC <- biomaRt::getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
       				            filters = c(filters = "entrezgene_id"),
       			         	    values = gene_list,
					    uniqueRows = TRUE,
			                    useCache = FALSE,
			                    mart = ensembl)

		#removing duplicated rows
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]
	
		# replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
		ID_to_HGNC[ID_to_HGNC == ""] <- NA
		ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC),]

		#add HGNC to dataframe by joining them
		dge_data <- dplyr::left_join(dge_data, ID_to_HGNC, 
       			                     by = c("GeneID" = 'entrezgene_id'),
                            		     copy = FALSE)

		#replace gene_id column with new identifier
		dge_data <- dge_data[, -1]
		dge_data <- dge_data[, c(7,1:6)]



	
	} else if (grepl("^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|XM|XP|XR|YP|ZP)_\\d+)|(NZ\\_[A-Z]{2,4}\\d+))(\\.\\d+)?$", map_data == TRUE)) {

		#REFSEQ
		ensembl <- biomaRt::useMart("ensembl")
		ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

		datasets <- biomaRt::listDatasets(ensembl)
		attributes <- biomaRt::listAttributes(ensembl)
		filters <- biomaRt::listFilters(ensembl)

		# prepare query identifiers from the transcripts dataset
		gene_list <- dge_data[ ,1]

		ID_to_HGNC <- biomaRt::getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
       				            filters = c(filters = "entrezgene_id"),
       			         	    values = gene_list,
					    uniqueRows = TRUE,
			                    useCache = FALSE,
			                    mart = ensembl)

		#removing duplicated rows
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
		ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]
	
		# replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
		ID_to_HGNC[ID_to_HGNC == ""] <- NA
		ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC),]

		#add HGNC to dataframe by joining them
		dge_data <- dplyr::left_join(dge_data, ID_to_HGNC, 
       			                     by = c("GeneID" = 'entrezgene_id'),
                            		     copy = FALSE)

		#replace gene_id column with new identifier
		dge_data <- dge_data[, -1]
		dge_data <- dge_data[, c(7,1:6)]

	} else {

		stop("Invalid identifier")
	}

}


#data import
if (is.character(dge_data) == TRUE) {
 
	if (grepl(".csv", as.character(dge_data)) == TRUE) {
    
		dge_data <- read.csv(dge_data) 

	} else if (grepl(".tabular", as.character(dge_data)) == TRUE) {

		dge_data <- read.table(dge_data)	
	
	} else {

		dge_data <- readxl::read_excel(dge_data, sheet = 1)

	}

}

if (type == "deseq2") {

	if (mapping == TRUE) {

		dge_data <- map(dge_data)

	}


	#data wrangling for visualization
	dge_data <- na.omit(dge_data)
	dge_data <- tibble::remove_rownames(dge_data)
	colnames(dge_data)[1] <- "genes"
	dge_data <- tibble::column_to_rownames(dge_data, var = "genes")
	colnames(dge_data)[1] <- "baseMean"
	colnames(dge_data)[2] <- "log2FoldChange"
	colnames(dge_data)[6] <- "padj"

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
