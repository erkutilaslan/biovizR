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
#' @param tech_rep Number of technical replicates. Necessary for accurate statistical calculation
#' @param test Default TRUE. Enable or disable statistical analysis.
#' @param stat Default t-test. Select statistical testing method.
#' @return A bar plot of qPCR results.
#' @import tidyverse
#' @export

#qpcr_data <- read.csv("/mnt/c/Users/Erkut Ilaslan/Desktop/FOXM1/RIP results/RIP_qPCR.csv")
#group1 <- "RIP NC"
#group2 <- "RIP P1"
#ref1 <- "Rluc"
#ref2 <- "Fluc"
#goi <- "FOXM1"
#tech_rep <- "4"

#qpcr_data <- read.csv("~/rip_pum1.csv")
#group1 <- "RIP NC"
#group2 <- "RIP P1"
#ref1 <- "Fluc"
#ref2 <- "Rluc"
#goi <- "FOXM1"
#tech_rep <- 3
#test <- TRUE

#qpcr_data <- read.csv("~/Cq siFOXM1 siPUM1.csv")
#type <- "Bio-rad"
#group1 <- "siCTRL 1"
#group2 <- "siFOXM1 1"
#ref1 <- "GARS1"
#ref2 <- ""
#goi <- "FOXM1"
#tech_rep <- 4
#test <- FALSE

#barplot_qpcr("~/Cq siFOXM1 siPUM1.csv",
#	     group1 = "siCTRL 1",
#	     group2 = "siFOXM1 1",
#	     ref1 = "GARS1",
#	     ref2 = "DTD1",
#	     goi = "FOXM1",
#	     tech_rep = 4,
#             test = FALSE)

barplot_qpcr <- function(qpcr_data,
                         type = "biorad",
			 group1 = "",
			 group2 = "",
                         ref1 = "",
			 ref2 = "",
                         goi = "",
			 tech_rep = 4,
			 test = TRUE,
                         stat = "t-test") {

#data import
qpcr_data <- read.csv(qpcr_data)

#data wrangling
if (type == "biorad") {

  if (test == TRUE) {

    qpcr_data <- qpcr_data[, c(3, 5, 6, 7)]

  } else {

    qpcr_data <- qpcr_data[, c(3, 5, 7)]

  }

  qpcr_data$Sample[qpcr_data$Sample == ""] <- NA
  qpcr_data$Target[qpcr_data$Target == "Target"] <- NA
  qpcr_data <- na.omit(qpcr_data)
  qpcr_data$Cq <- as.numeric(as.character(qpcr_data$Cq))

  if (test == TRUE) {

    qpcr_data <- dplyr::group_by(qpcr_data, Target, Sample, Biological.Set.Name)

  } else {

    qpcr_data <- dplyr::group_by(qpcr_data, Target, Sample)

  }

  qpcr_data <- dplyr::mutate(qpcr_data, Cq_mean = mean(Cq))
  qpcr_data <- dplyr::ungroup(qpcr_data)

  if (test == TRUE) {

    qpcr_data <- qpcr_data[ ,-4]

  } else {

    qpcr_data <- qpcr_data[ ,-3]

  }

  qpcr_data <- dplyr::distinct(qpcr_data)

#pivot_wider for analysis
  qpcr_data <- tidyr::pivot_wider(qpcr_data,
                                  names_from = Target, values_from = Cq_mean)

}

#analysis
#for each sample average the ct values of the 3 ref genes --> ct[ref] = mean( ct[ref1], ct[ref[2], ct[ref3] )
#for each sample, calculate the dct as the difference dct = ct[ref] - ct[goi]

#annotating ref and goi names
colnames(qpcr_data)[grep(ref1, colnames(qpcr_data))] <- "ref1"

if (ref2 != "") {

  colnames(qpcr_data)[grep(ref2, colnames(qpcr_data))] <- "ref2"
}

colnames(qpcr_data)[grep(goi, colnames(qpcr_data))] <- "goi"

#avg of refs for multiple refs
if (ref2 != "") {

  qpcr_data <- dplyr::mutate(qpcr_data, avg_ref = (ref1 + ref2) / 2)
#dCt calculation for avg of refs
  qpcr_data <- dplyr::mutate(qpcr_data, dCt = (goi - avg_ref))

} else {

#dCt calculation for only 1 ref

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
qpcr_data <- dplyr::left_join(qpcr_data, ref_exp)
qpcr_data <- dplyr::left_join(qpcr_data, target_exp)
qpcr_data <- dplyr::mutate(qpcr_data, avg_exp = dplyr::coalesce(target_avg_exp, ref_avg_exp))

#converting raw expression to percentage for better visualization
qpcr_data <- dplyr::arrange(qpcr_data, dplyr::desc(qpcr_data$ref_avg_exp))
qpcr_data <- dplyr::mutate(qpcr_data, percent_exp = qpcr_data$avg_exp/qpcr_data$ref_avg_exp[1]*100)

if (test == TRUE) {

#duplicating rows for accurate p-value calculation
  idx <- rep(1:nrow(target_exp), tech_rep)
  target_exp <- target_exp[idx, ]

  idx2 <- rep(1:nrow(ref_exp), tech_rep)
  ref_exp <- ref_exp[idx2, ]

#calculating p-value
  stats <- t.test(target_exp$expression, ref_exp$expression, alternative = "two.sided")

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

}

#calculating sd
target_sd <- sd(target_exp$expression)
ref_sd <- sd(ref_exp$expression)

#converting raw sd into percentage
target_sd <- target_sd/qpcr_data$target_avg_exp[4]*100
ref_sd <- ref_sd/qpcr_data$ref_avg_exp[1]*100

#distinct only 1 column
qpcr_data <- qpcr_data[!duplicated(qpcr_data$percent_exp), ]

#removing na to only visualize sample of interest
qpcr_data <- qpcr_data[!is.na(qpcr_data$percent_exp), ]

#visualization
if (test == TRUE) {

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
					 y.position = qpcr_data$percent_exp[2] + ref_sd + 10,
					 bracket.size = 0.6, label.size = 6)

} else {

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
     					 width = 0.1))

}

plot(final_plot)
return(final_plot)

}

