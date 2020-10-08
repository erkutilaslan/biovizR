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
actb <- read.csv("~/ACTB.csv")
gapdh <- read.csv("~/GAPDH.csv")
nanos1 <- read.csv("~/NANOS1.csv")

#data wrangling
#Detector = Gene 1, Gene2, ...
#Platename = name of the Plate ran for analysis
#Sample is sample the same

#merging all files
qpcr_data <- dplyr::bind_rows(actb, gapdh, nanos1)
#naming columns
colnames(qpcr_data)[8] <- "Ct"
colnames(qpcr_data)[4] <- "Detector"
colnames(qpcr_data)[3] <- "Platename"
qpcr_data <- qpcr_data[, c(-1, -2, -7, -9:-16)]
qpcr_data <- na.omit(qpcr_data)

#setting parameters
ctrlsamples <- c("- control 24", "- control 48", "- control 72")
hkgenes <- c("GAPDH", "Actin")

#analysis
qpcr_results <- ddCtExpression(qpcr_data,
                            calibrationSample = ctrlsamples,
                            housekeepingGenes = hkgenes)

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
