#'
#' This function visualizes qPCR results as bar plot.
#'
#' @param qpcr_data Input qPCR data to visualize.
#' @param type Default Bio-rad. Specify model of thermocycler.
#' @param group1 Control group.
#' @param group2 Target group.
#' @param ref1 Reference gene one.
#' @param ref2 Reference gene two.
#' @param goi Gene of interest.
#' @param stat Default t-test. Select statistical testing method.
#' @return A bar plot of qPCR results.
#' @import tidyverse
#' @export


barplot_qpcr <- function(qpcr_data,
                         type = "biorad",
			 group1 = "",
			 group2 = "",
                         ref1 = "",
			 ref2 = "",
                         goi = "",
                         stat = "t-test") {

#data import
qpcr_data <- read.csv(qpcr_data)

#data wrangling
if (type == "biorad") {

qpcr_data <- qpcr_data[, c(3, 5, 6, 7)]
qpcr_data$Sample[qpcr_data$Sample == ""] <- NA
qpcr_data$Target[qpcr_data$Target == "Target"] <- NA
qpcr_data <- na.omit(qpcr_data)
qpcr_data <- dplyr::distinct(qpcr_data)
qpcr_data$Cq <- as.numeric(as.character(qpcr_data$Cq))
qpcr_data <- dplyr::group_by(qpcr_data, Sample, Biological.Set.Name)
qpcr_data <- dplyr::mutate(qpcr_data, Cq_mean = mean(Cq))

#pivot_wider for analysis
qpcr_data <- tidyr::pivot_wider(qpcr_data, names_from = Target, values_from = Cq)

}

#analysis
#for each sample average the ct values of the 3 ref genes --> ct[ref] = mean( ct[ref1], ct[ref[2], ct[ref3] )
#for each sample, calculate the dct as the difference dct = ct[ref] - ct[goi]

#annotating ref and goi names
colnames(qpcr_data)[grep(ref1, colnames(qpcr_data))] <- "ref1"
colnames(qpcr_data)[grep(ref2, colnames(qpcr_data))] <- "ref2"
colnames(qpcr_data)[grep(goi, colnames(qpcr_data))] <- "goi"

#avg of refs for multiple refs
if (ref2 != "") {

qpcr_data <- dplyr::mutate(qpcr_data, avg_ref = (ref1 + ref2) / 2)

#dCt calculation for avg of refs
qpcr_data <- dplyr::mutate(qpcr_data, dCt = (goi - avg_ref))

}

#dCt calculation for only 1 ref
if (ref2 == "") {
qpcr_data <- dplyr::mutate(qpcr_data, dCt = (goi - ref1))
}

#avg.dCt of target gene value in control biological replicates
ref_ct <- qpcr_data[grep(group1, qpcr_data$Sample), ]
qpcr_data <- dplyr::mutate(qpcr_data, avg_ref_dCt = mean(ref_ct$goi))

#ddCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, ddCt = dCt - avg_ref_dCt)

#2^-ddCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, expression = 2^-ddCt)

#mean the expression in biological replicates
ref_exp <- qpcr_data[grep(group1, qpcr_data$Sample), ]
target_exp <- qpcr_data[grep(group2, qpcr_data$Sample), ]
ref_exp <- dplyr::mutate(ref_exp, ref_avg_exp = mean(expression))
target_exp <- dplyr::mutate(target_exp, target_avg_exp = mean(expression))
qpcr_data <- plyr::join_all(list(qpcr_data, ref_exp, target_exp))
qpcr_data <- dplyr::mutate(qpcr_data, avg_exp = dplyr::coalesce(target_avg_exp, ref_avg_exp))

#converting raw expression to percentage for better visualization
qpcr_data <- dplyr::arrange(qpcr_data, dplyr::desc(qpcr_data$avg_exp))
qpcr_data <- dplyr::mutate(qpcr_data, percent_exp = qpcr_data$avg_exp/qpcr_data$avg_exp[1]*100)

#calculating p-value
stats <- t.test(target_exp$goi, ref_exp$goi)

#converting pvalues to *
if (stats$p.value > 0.05) {
	pvalue <- "ns"
}

if (stats$p.value < 0.05) {
	pvalue <- "*"
}

if (stats$p.value < 0.01) {
	pvalue <- "**"
}

if (stats$p.value < 0.001) {
	pvalue <- "***"
}

#we need to change group names to represent the data automatically
#this is the supported df layout for ggpubr::stat_pvalue_manuel()
qpcr_data2 <- tibble::tribble(~group1, ~group2, ~pvalue,
                              group1, group2, pvalue)

#calculating sd
target_sd <- sd(target_exp$expression)
ref_sd <- sd(ref_exp$expression)

#qpcr_data <- qpcr_data[c(-2,-3,-5,-6), ]

#converting raw sd into percentage
target_sd <- target_sd/qpcr_data$avg_exp[1]*100
ref_sd <- ref_sd/qpcr_data$avg_exp[1]*100

#distinct only 1 column
qpcr_data <- qpcr_data[!duplicated(qpcr_data$percent_exp), ]

#visualization
final_plot <- ggpubr::ggbarplot(qpcr_data,
				x = "Sample",
				y = "percent_exp",
				fill = "808080",
				xlab = "Sample",
				ylab = "Relative mRNA level",
				size = 0.5,
				palette = "npg",
				lab.size = 5,
				lab.vjust = 0.5,
				lab.hjust = 1.2,
				sort.by.groups = FALSE,
				ggtheme = ggpubr::theme_pubr(base_size = 14)) +
              ggplot2::geom_errorbar(ggplot2::aes(x = group2,
	        			 ymin = percent_exp[2] - target_sd,
	                		 ymax = percent_exp[2] + target_sd,
		                	 width = 0.1)) +
              ggplot2::geom_errorbar(ggplot2::aes(x = group1,
	    				 ymin = percent_exp[1] - ref_sd,
     					 ymax = percent_exp[1] + ref_sd,
     					 width = 0.1)) +
              ggpubr::stat_pvalue_manual(qpcr_data2, label = "pvalue",
					 y.position = qpcr_data$percent_exp[1] + ref_sd + 10,
					 bracket.size = 0.6)

plot(final_plot)
return(final_plot)

}
