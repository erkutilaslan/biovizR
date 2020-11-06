#'
#' This function visualizes qPCR results as bar plot.
#'
#' @param qpcr_data Input qPCR data to visualize.
#' @param type 
#' @param ref 
#' @param goi
#' @param stat Default t-test.
#' @return A bar plot of qPCR results.
#' @import tidyverse
#' @export


barplot_qpcr <- function(qpcr_data,
                         type = "biorad",
                         ref = "",
                         goi = "",
                         stat = "t-test") {

#data import
qpcr_data <- read.csv(qpcr_data)
qpcr_data <- read.csv("~/biovizR_data/qpcr_data.csv")

#data wrangling
qpcr_data <- qpcr_data[, c(4, 6, 9)]
qpcr_data$Sample[qpcr_data$Sample == ""] <- NA
qpcr_data$Target[qpcr_data$Target == "Target"] <- NA
qpcr_data <- na.omit(qpcr_data)
qpcr_data <- dplyr::distinct(qpcr_data)
qpcr_data$Cq.Mean <- as.numeric(as.character(qpcr_data$Cq.Mean))

#pivot_wider for analysis
qpcr_data <- tidyr::pivot_wider(qpcr_data, names_from = Target, values_from = Cq.Mean)

#remove so that we can write the analysis script properly will be removed later on
qpcr_data <- qpcr_data[c(-2, -5, -8), ]

#analysis
#for each sample average the ct values of the 3 ref genes --> ct[ref] = mean( ct[ref1], ct[ref[2], ct[ref3] )
#for each sample, calculate the dct as the difference dct = ct[ref] - ct[goi]

#here avg of refs for multiple refs

#dCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, dCt = (NOS1 - Actin))

#avg.dCt of target gene value in control biological replicates
ref_ct <- qpcr_data[grep("control", qpcr_data$Sample), ]
qpcr_data <- dplyr::mutate(qpcr_data, avg_ref_dCt = mean(ref_ct$NOS1))

#ddCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, ddCt = dCt - avg_ref_dCt)

#2^-ddCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, expression = 2^-ddCt)

#mean the expression in biological replicates
ref_exp <- qpcr_data[grep("control", qpcr_data$Sample), ]
target_exp <- qpcr_data[grep("N1-1", qpcr_data$Sample), ]
ref_exp <- dplyr::mutate(ref_exp, ref_avg_exp = mean(expression))
target_exp <- dplyr::mutate(target_exp, target_avg_exp = mean(expression))
qpcr_data <- plyr::join_all(list(qpcr_data, ref_exp, target_exp))
qpcr_data <- dplyr::mutate(qpcr_data, avg_exp = dplyr::coalesce(target_avg_exp, ref_avg_exp))

#here try to generate qpcr_data2 in a way that t_test function from rstatix can process and calculate
#try tibble::tribble or any other method to transform the df in rstatix compatible format

#converting raw expression to percentage for better visualization
qpcr_data <- dplyr::arrange(qpcr_data, dplyr::desc(qpcr_data$avg_exp))
qpcr_data <- dplyr::mutate(qpcr_data, percent_exp = qpcr_data$avg_exp/qpcr_data$avg_exp[1]*100)

#calculating p-value
stats <- t.test(target_exp$NOS1, ref_exp$NOS1)
#if statement for if stats$pvalue < 0.05 = *, if <0.01 **, if pvalue < 0.001 ***
pvalue <- stats$p.value
qpcr_data2 <- tibble::tribble(~group1, ~group2, ~pvalue,
                              "N1-1 24", "- control 24", pvalue)

#calculating sd
target_sd <- sd(target_exp$expression)
ref_sd <- sd(ref_exp$expression)

qpcr_data <- qpcr_data[c(-2,-3,-5,-6), ]

#converting raw sd into percentage
target_sd <- target_sd/qpcr_data$avg_exp[1]*100
ref_sd <- ref_sd/qpcr_data$avg_exp[1]*100

#visualization
final_plot <- ggpubr::ggbarplot(qpcr_data,
				x = "Sample",
				y = "percent_exp",
				fill = "Sample",
				xlab = "Sample",
				ylab = "Relative mRNA level",
				size = 0.5,
				palette = "npg",
				lab.size = 5,
				lab.vjust = 0.5,
				lab.hjust = 1.2,
				sort.by.groups = FALSE,
				ggtheme = ggpubr::theme_pubr(base_size = 14)) +
              ggplot2::geom_errorbar(ggplot2::aes(x = "N1-1 24",
	        			 ymin = percent_exp[2] - target_sd,
	                		 ymax = percent_exp[2] + target_sd,
		                	 width = 0.1)) +
              ggplot2::geom_errorbar(ggplot2::aes(x = "- control 24",
	    				 ymin = percent_exp[1] - ref_sd,
     					 ymax = percent_exp[1] + ref_sd,
     					 width = 0.1)) +
              ggpubr::stat_pvalue_manual(qpcr_data2, label = "pvalue",
					 y.position = qpcr_data$percent_exp[1] + ref_sd + 10,
					 bracket.size = 0.6)

plot(final_plot)
return(final_plot)

}
