#' Dot plot for IP-MS experiment visualization.
#'
#' This function visualizes IP-MS results with 2 axis
#' for number of validated peptides in 2 samples.
#'
#' @param ip_data Path to the input file
#' @return A dot plot of the ip_data
#' @import tidyverse
#' @export

dotplot_ipms <- function(ip_data) {

#data import
ip_1 <- readxl::read_xlsx(ip_data,
                  sheet = 3)
ip_2 <- readxl::read_xlsx(ip_data,
                  sheet = 4)
ip_3 <- readxl::read_xlsx(ip_data,
                  sheet = 5)

ctrl_1 <- readxl::read_xlsx(ip_data,
                  sheet = 6)
ctrl_2 <- readxl::read_xlsx(ip_data,
                  sheet = 7)
ctrl_3 <- readxl::read_xlsx(ip_data,
                  sheet = 8)

#data cleanup and process

#removing unnecesary columns
ip_1 <- ip_1[, c(4, 8)]
ip_2 <- ip_2[, c(4, 8)]
ip_3 <- ip_3[, c(4, 8)]
ctrl_1 <- ctrl_1[, c(4, 8)]
ctrl_2 <- ctrl_2[, c(4, 8)]
ctrl_3 <- ctrl_3[, c(4, 8)]

#renaming columns for join
colnames(ip_1) <- c("HGNC", "ip_1")
colnames(ip_2) <- c("HGNC", "ip_2")
colnames(ip_3) <- c("HGNC", "ip_3")
colnames(ctrl_1) <- c("HGNC", "ctrl_1")
colnames(ctrl_2) <- c("HGNC", "ctrl_2")
colnames(ctrl_3) <- c("HGNC", "ctrl_3")

#joining DFs into one dataframe
ipms <- dplyr::full_join(ip_1, ip_2, "HGNC")
ipms <- dplyr::full_join(ipms, ip_3, "HGNC")
ipms <- dplyr::full_join(ipms, ctrl_1, "HGNC")
ipms <- dplyr::full_join(ipms, ctrl_2, "HGNC")
ipms <- dplyr::full_join(ipms, ctrl_3, "HGNC")

#dropna to remove NA from only one column
ipms <- tidyr::drop_na(ipms, HGNC)

#replace NA in peptide numbers with 0
ipms <- tidyr::replace_na(ipms, list(ip_1 = 0,
                              ip_2 = 0,
                              ip_3 = 0,
                              ctrl_1 = 0,
                              ctrl_2 = 0,
                              ctrl_3 = 0))

#taking the average of peptide numbers using mutate
ipms <- dplyr::mutate(ipms, "avg.ip" = (ip_1 + ip_2 + ip_3) / 3)
ipms <- dplyr::mutate(ipms, "avg.ctrl" =  (ctrl_1 + ctrl_2 + ctrl_3) / 3)

#generating ma-plot

ip_plot <- ggplot2::ggplot(ipms, ggplot2::aes(avg.ip, avg.ctrl), label = HGNC)
final_plot <- ip_plot + ggplot2::geom_point() +
    ggplot2::geom_text(aes(label = ifelse(avg.ctrl < 2,
                                 as.character(HGNC),
                                 "")),
                  hjust = 0, vjust = 0)

plot(final_plot)
return(final_plot)
}

