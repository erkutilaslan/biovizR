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


barplot_qpcr <- function(qpcr_data,
                         type = "biorad",
                         ref = "",
                         target = "",
                         stat = "t-test") {

#data import
qpcr_data <- read.csv(qpcr_data)
qpcr_data <- read.csv("~/qpcr_data.csv")

#data wrangling
qpcr_data <- qpcr_data[, c(4, 6, 9)]
qpcr_data$Sample[qpcr_data$Sample == ""] <- NA
qpcr_data$Target[qpcr_data$Target == "Target"] <- NA
qpcr_data <- na.omit(qpcr_data)
qpcr_data <- dplyr::distinct(qpcr_data)
qpcr_data$Cq.Mean <- as.numeric(as.character(qpcr_data$Cq.Mean))

#pivot_wider for analysis
qpcr_data <- tidyr::pivot_wider(qpcr_data, names_from = Target, values_from = Cq.Mean )

#remove so that we can write the analysis script properly will be removed later on
qpcr_data <- qpcr_data[c(-2, -5, -8), ]

#analysis
#for each sample average the ct values of the 3 ref genes --> ct[ref] = mean( ct[ref1], ct[ref[2], ct[ref3] )
# for each sample, calculate the dct as the difference dct = ct[ref] - ct[goi]

#dCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, dCt = (NOS1 - Actin))

#avg.dCt of target gene value in control biological replicates
control_ct <- qpcr_data[grep("control", qpcr_data$Sample), ]
qpcr_data <- dplyr::mutate(qpcr_data, avg.dCt = mean(control_ct$NOS1))

#ddCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, ddCt = dCt - avg.dCt)

#2^-ddCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, expression = 2^-ddCt)

#mean the expression in biological replicates
control_exp <- qpcr_data[grep("control", qpcr_data$Sample), ]
target_exp <- qpcr_data[grep("N1-1", qpcr_data$Sample), ]
control_exp <- dplyr::mutate(control_exp, control_avg.exp = mean(expression))
target_exp <- dplyr::mutate(target_exp, target_avg.exp = mean(expression))
qpcr_data <- plyr::join_all(list(qpcr_data, control_exp, target_exp))
qpcr_data <- dplyr::mutate(qpcr_data, avg.exp = dplyr::coalesce(qpcr_data$target_avg.exp, control_avg.exp))

#statistical analysis
stats <- t.test(target_exp$NOS1, control_exp$NOS1)
stats$p.value

#visualization
final_plot <- ggpubr::ggbarplot(test2,
          x = "Sample",
          y = "avg.exp",
#          fill = "darkgray",
          xlab = "Sample",
          ylab = "Relative mRNA level",
          size = 0.5,
          palette = "jco",
          title = "RT-qPCR",
          lab.size = 5,
          lab.vjust = 0.5,
          lab.hjust = 1.2,
          sort.by.groups = FALSE,
          ggtheme = ggpubr::theme_pubr(base_size = 18))

plot(final_plot)
return(final_plot)

}
