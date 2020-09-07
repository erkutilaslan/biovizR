#' Dot plot for IP-MS experiment visualization.
#'
#' This function visualizes IP-MS results with 2 axis
#' for number of validated peptides in NANOS-IP vs pCMV6-IP
#' when i run the script manually it generates proper plot
#' but when i run it as a function, it proceeds without error
#' and creates pdf file but the file doesnt open.
#'
#' @param ip_data Path to the input file
#' @return A dot plot of the ip_data
#' @export

dotplot_ipms <- function(ip_data) {

#data import
ip_1 <- read_xlsx(ip_data,
                  sheet = 3)
ip_2 <- read_xlsx(ip_data,
                  sheet = 4)
ip_3 <- read_xlsx(ip_data,
                  sheet = 5)

ctrl_1 <- read_xlsx(ip_data,
                  sheet = 6)
ctrl_2 <- read_xlsx(ip_data,
                  sheet = 7)
ctrl_3 <- read_xlsx(ip_data,
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
ipms <- full_join(ip_1, ip_2, "HGNC")
ipms <- full_join(ipms, ip_3, "HGNC")
ipms <- full_join(ipms, ctrl_1, "HGNC")
ipms <- full_join(ipms, ctrl_2, "HGNC")
ipms <- full_join(ipms, ctrl_3, "HGNC")

#dropna to remove NA from only one column
ipms <- ipms %>% drop_na(HGNC)

#replace NA in peptide numbers with 0
ipms <- replace_na(ipms, list(ip_1 = 0,
                              ip_2 = 0,
                              ip_3 = 0,
                              ctrl_1 = 0,
                              ctrl_2 = 0,
                              ctrl_3 = 0))

#taking the average of peptide numbers using mutate
ipms <- ipms %>% mutate("avg.ip" = (ip_1 + ip_2 + ip_3) / 3)
ipms <- ipms %>% mutate("avg.ctrl" =  (ctrl_1 + ctrl_2 + ctrl_3) / 3)

#generating ma-plot

ip_plot <- ggplot(ipms, aes(avg.ip, avg.ctrl), label = HGNC)
final_plot <- ip_plot + geom_point() +
    geom_text(aes(label = ifelse(avg.ctrl < 2,
                                 as.character(HGNC),
                                 "")),
                  hjust = 0, vjust = 0)

plot(final_plot)
return(final_plot)
}

