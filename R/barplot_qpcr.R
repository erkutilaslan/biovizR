#' Analyze and visualize qPCR results with barplot.
#'
#' This function visualizes qPCR results as bar plot.
#'
#' @param qpcr_data Input qPCR data to visualize.
#' @param type 
#' @param ref 
#' @param target
#' @param stat Default t-test.
#' @return A bar plot of qPCR results.
#' @import tidyverse
#' @export

library(ddCt)
barplot_qpcr <- function(qpcr_data,
                         type = "biorad",
                         ref = "",
                         target = "",
                         stat = "t-test") {

#data import
qpcr_data <- readxl::read_xlsx(qpcr_data, sheet = 1)

actb <- read.csv("~/biovizR_data/ACTB.csv")
gapdh <- read.csv("~/biovizR_data/GAPDH.csv")
nanos1 <- read.csv("~/biovizR_data/NANOS1.csv")

#data wrangling

#merging all files
qpcr_data <- dplyr::bind_rows(actb, gapdh, nanos1)

#naming columns
qpcr_data <- qpcr_data[, c(-1, -2, -7, -9:-16)]
qpcr_data <- na.omit(qpcr_data)

#analysis
#for each sample average the ct values of the 3 ref genes --> ct[ref] = mean( ct[ref1], ct[ref[2], ct[ref3] )
# for each sample, calculate the dct as the difference dct = ct[ref] - ct[goi]


#visualization
final_plot <- ggpubr::ggbarplot(qpcr_results,
          x = "Samples",
          y = "Ct mean",
          fill = "darkgray",
          xlab = "Samples",
          ylab = "Relative mRNA level",
          size = 0.5,
          palette = "jco",            # jco journal color palett. see ?ggpar
          title = "RT-qPCR",
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          sort.by.groups = FALSE,     # Don't sort inside each group
          ggtheme = ggpubr::theme_pubr(base_size = 18))

plot(final_plot)
return(final_plot)

}
