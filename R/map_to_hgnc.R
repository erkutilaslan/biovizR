#' Auto mapping to HGNC of a variety of gene identifiers.
#'
#' This function automatically detects a variety of gene identifiers and converts them to HGNC using biomaRt.
#'
#' @param map_data Path to the input file.
#' @return A table containing HGNC gene IDs.
#' @export

map_data <- "~/Work/PhD/Transcriptomics/N1 N3 OE RNA-Seq/NANOS1_1_counts.tabular"
dat <- read.table(map_data)

# mapping function
map_to_hgnc <- function(map_data) {

  # data import
  if (is.character(map_data) == TRUE) {
    if (grepl(".csv", as.character(map_data)) == TRUE) {
      map_data <- read.csv(map_data, header = TRUE)
    } else {
      map_data <- readxl::read_excel(map_data, sheet = 1)
    }
  }

  if (grepl("^\\d+$", map_data[1, 1]) == TRUE) {
    colnames(map_data)[1] <- "GeneID"
    ensembl <- biomaRt::useMart("ensembl")
    ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

    datasets <- biomaRt::listDatasets(ensembl)
    attributes <- biomaRt::listAttributes(ensembl)
    filters <- biomaRt::listFilters(ensembl)

    # prepare query identifiers from the transcripts dataset
    gene_list <- map_data[, 1]

    ID_to_HGNC <- biomaRt::getBM(
      attributes = c("entrezgene_id", "hgnc_symbol"),
      filters = c(filters = "entrezgene_id"),
      values = gene_list,
      uniqueRows = TRUE,
      useCache = FALSE,
      mart = ensembl
    )

    # removing duplicated rows
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]

    # replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
    ID_to_HGNC[ID_to_HGNC == ""] <- NA
    ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC), ]

    # add HGNC to dataframe by joining them
    map_data <- dplyr::left_join(map_data, ID_to_HGNC,
      by = c("GeneID" = "entrezgene_id"),
      copy = FALSE
    )

    # replace gene_id column with new identifier
    map_data <- map_data[, -1]
    map_data <- map_data[, c(ncol(map_data), 1:ncol(map_data))]
    map_data <- map_data[, -ncol(map_data)]

  } else if (grepl("^(\\w+\\d+(\\.\\d+)?)|(NP_\\d+)$", map_data[1, 1]) == TRUE) {

    # NCBI Protein
    ensembl <- biomaRt::useMart("ensembl")
    ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

    datasets <- biomaRt::listDatasets(ensembl)
    attributes <- biomaRt::listAttributes(ensembl)
    filters <- biomaRt::listFilters(ensembl)

    # prepare query identifiers from the transcripts dataset
    gene_list <- map_data[, 1]

    ID_to_HGNC <- biomaRt::getBM(
      attributes = c("entrezgene_id", "hgnc_symbol"),
      filters = c(filters = "entrezgene_id"),
      values = gene_list,
      uniqueRows = TRUE,
      useCache = FALSE,
      mart = ensembl
    )

    # removing duplicated rows
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]

    # replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
    ID_to_HGNC[ID_to_HGNC == ""] <- NA
    ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC), ]

    # add HGNC to dataframe by joining them
    map_data <- dplyr::left_join(map_data, ID_to_HGNC,
      by = c("GeneID" = "entrezgene_id"),
      copy = FALSE
    )

    # replace gene_id column with new identifier
    map_data <- map_data[, -1]
    map_data <- map_data[, c(7, 1:6)]
  } else if (grepl("^RF\\d{5}$", map_data[1, 1]) == TRUE) {

    # RFAM
    ensembl <- biomaRt::useMart("ensembl")
    ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

    datasets <- biomaRt::listDatasets(ensembl)
    attributes <- biomaRt::listAttributes(ensembl)
    filters <- biomaRt::listFilters(ensembl)

    # prepare query identifiers from the transcripts dataset
    gene_list <- map_data[, 1]

    ID_to_HGNC <- biomaRt::getBM(
      attributes = c("entrezgene_id", "hgnc_symbol"),
      filters = c(filters = "entrezgene_id"),
      values = gene_list,
      uniqueRows = TRUE,
      useCache = FALSE,
      mart = ensembl
    )

    # removing duplicated rows
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]

    # replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
    ID_to_HGNC[ID_to_HGNC == ""] <- NA
    ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC), ]

    # add HGNC to dataframe by joining them
    map_data <- dplyr::left_join(map_data, ID_to_HGNC,
      by = c("GeneID" = "entrezgene_id"),
      copy = FALSE
    )

    # replace gene_id column with new identifier
    map_data <- map_data[, -1]
    map_data <- map_data[, c(7, 1:6)]
  } else if (grepl("^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\\.\\d+)?$", map_data == TRUE)) {

    # UNIPROT
    ensembl <- biomaRt::useMart("ensembl")
    ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

    datasets <- biomaRt::listDatasets(ensembl)
    attributes <- biomaRt::listAttributes(ensembl)
    filters <- biomaRt::listFilters(ensembl)

    # prepare query identifiers from the transcripts dataset
    gene_list <- map_data[, 1]

    ID_to_HGNC <- biomaRt::getBM(
      attributes = c("entrezgene_id", "hgnc_symbol"),
      filters = c(filters = "entrezgene_id"),
      values = gene_list,
      uniqueRows = TRUE,
      useCache = FALSE,
      mart = ensembl
    )

    # removing duplicated rows
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]

    # replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
    ID_to_HGNC[ID_to_HGNC == ""] <- NA
    ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC), ]

    # add HGNC to dataframe by joining them
    map_data <- dplyr::left_join(map_data, ID_to_HGNC,
      by = c("GeneID" = "entrezgene_id"),
      copy = FALSE
    )

    # replace gene_id column with new identifier
    map_data <- map_data[, -1]
    map_data <- map_data[, c(7, 1:6)]
  } else if (grepl("^((ENS[FPTG]\\d{11}(\\.\\d+)?)|(FB\\w{2}\\d{7})|(Y[A-Z]{2}\\d{3}[a-zA-Z](\\-[A-Z])?)|([A-Z_a-z0-9]+(\\.)?(t)?(\\d+)?([a-z])?))$", map_data == TRUE)) {

    # ENSEMBL
    ensembl <- biomaRt::useMart("ensembl")
    ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

    datasets <- biomaRt::listDatasets(ensembl)
    attributes <- biomaRt::listAttributes(ensembl)
    filters <- biomaRt::listFilters(ensembl)

    # prepare query identifiers from the transcripts dataset
    gene_list <- map_data[, 1]

    ID_to_HGNC <- biomaRt::getBM(
      attributes = c("entrezgene_id", "hgnc_symbol"),
      filters = c(filters = "entrezgene_id"),
      values = gene_list,
      uniqueRows = TRUE,
      useCache = FALSE,
      mart = ensembl
    )

    # removing duplicated rows
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]

    # replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
    ID_to_HGNC[ID_to_HGNC == ""] <- NA
    ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC), ]

    # add HGNC to dataframe by joining them
    map_data <- dplyr::left_join(map_data, ID_to_HGNC,
      by = c("GeneID" = "entrezgene_id"),
      copy = FALSE
    )

    # replace gene_id column with new identifier
    map_data <- map_data[, -1]
    map_data <- map_data[, c(7, 1:6)]
  } else if (grepl("^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|XM|XP|XR|YP|ZP)_\\d+)|(NZ\\_[A-Z]{2,4}\\d+))(\\.\\d+)?$", map_data == TRUE)) {

    # REFSEQ
    ensembl <- biomaRt::useMart("ensembl")
    ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

    datasets <- biomaRt::listDatasets(ensembl)
    attributes <- biomaRt::listAttributes(ensembl)
    filters <- biomaRt::listFilters(ensembl)

    # prepare query identifiers from the transcripts dataset
    gene_list <- map_data[, 1]

    ID_to_HGNC <- biomaRt::getBM(
      attributes = c("entrezgene_id", "hgnc_symbol"),
      filters = c(filters = "entrezgene_id"),
      values = gene_list,
      uniqueRows = TRUE,
      useCache = FALSE,
      mart = ensembl
    )

    # removing duplicated rows
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$entrezgene_id), ]
    ID_to_HGNC <- ID_to_HGNC[!duplicated(ID_to_HGNC$hgnc_symbol), ]

    # replace un-mapped empty cells with NA's, then drop all incomplete annotation cases.
    ID_to_HGNC[ID_to_HGNC == ""] <- NA
    ID_to_HGNC <- ID_to_HGNC[complete.cases(ID_to_HGNC), ]

    # add HGNC to dataframe by joining them
    map_data <- dplyr::left_join(map_data, ID_to_HGNC,
      by = c("GeneID" = "entrezgene_id"),
      copy = FALSE
    )

    # replace gene_id column with new identifier
    map_data <- map_data[, -1]
    map_data <- map_datmap_data[, c(7, 1:6)]
  } else {
    stop("Invalid identifier")
  }
}
