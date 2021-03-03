#'
#' This function visualizes qPCR results as bar plot.
#'
#' @param qpcr_data Input qPCR data to visualize.
#' @param type Default biorad. Specify model of thermocycler.
#' @param group1 Control group.
#' @param group2 Target group.
#' @param group3 Target group.
#' @param group4 Target group.
#' @param group5 Target group.
#' @param ref1 Reference gene one.
#' @param ref2 Reference gene two.
#' @param goi Gene of interest.
#' @param tech_rep Number of technical replicates. Optional if consideration of technical replicates for statistical testing is desired.
#' @param test Default TRUE. Enable or disable statistical analysis.
#' @param stat Default t-test. Select statistical testing method.
#' @return A bar plot of qPCR results.
#' @import tidyverse
#' @export

#qpcr_data <- read.csv("~/qpcr_foxm1.csv")
#group1 <- "siCTRL"
#group2 <- "siPUM1"
#group3 <- "siNANOS3"
#group4 <- "siPUM1/NANOS3"
#group5 <- "siFOXM1"
#ref1 <- "GARS1"
#ref2 <- ""
#goi <- "FOXM1"
#tech_rep <- "3"
#type <- "biorad"
#test <- FALSE

#qpcr_data <- read.csv("~/RIP_qPCR.csv")
#group1 <- "RIP NC"
#group2 <- "RIP P1"
#group3 <- "RIP N3"
#group4 <- ""
#group5 <- ""
#ref1 <- "Fluc"
#ref2 <- "Rluc"
#goi <- "FOXM1"
#tech_rep <- 3
#test <- FALSE
#type <- "biorad"
#stat <- "t.test"

#qpcr_data <- read.csv("~/Cq siFOXM1 siPUM1.csv")
#type <- "biorad"
#group1 <- "siCTRL 1"
#group2 <- "siFOXM1 1"
#ref1 <- "GARS1"
#ref2 <- ""
#goi <- "FOXM1"
#tech_rep <- 4
#test <- FALSE

#qpcr_data <- read.csv("~/multipe_test_qpcr.csv")
#type <- "biorad"
#group1 <- "5"
#group2 <- "6"
#group3 <- ""
#group4 <- ""
#group5 <- ""
#ref1 <- "GAPDH"
#ref2 <- ""
#goi <- "DAZL"
#tech_rep <- 3
#test <- FALSE

#barplot_qpcr("~/qpcr_foxm1.csv",
#             group1 = "siCTRL",
#             group2 = "siPUM1",
#	     group3 = "siNANOS3",
#	     group4 = "siPUM1/NANOS3",
#	     group5 = "siFOXM1",
#             ref1 = "GARS1",
#             goi = "FOXM1",
#	     tech_rep = 3,
#             test = FALSE)

#barplot_qpcr("~/multipe_test_qpcr.csv",
#             group1 = "5",
#             group2 = "6",
#	     group3 = "7",
#	     group4 = "8",
#	     group5 = "9",
#             ref1 = "GAPDH",
#             goi = "DAZL",
#             tech_rep = 3,
#             test = FALSE)


#barplot_qpcr("~/rip_pum1.csv",
 #            group1 = "RIP NC",
  #           group2 = "RIP P1",
   #          ref1 = "Rluc",
    #         ref2 = "Fluc",
     #        goi = "FOXM1",
      #       tech_rep = 3,
       #      test = FALSE)

barplot_qpcr <- function(qpcr_data,
                         type = "biorad",
			 group1 = "",
			 group2 = "",
			 group3 = "",
			 group4 = "",
			 group5 = "",
                         ref1 = "",
			 ref2 = "",
                         goi = "",
			 tech_rep = 3,
			 test = TRUE,
                         stat = "t-test") {

#data import
if (is.character(qpcr_data) == TRUE) {

  qpcr_data <- read.csv(qpcr_data)

}

#data wrangling
if (type == "biorad") {

  if (colnames(qpcr_data[1]) == "X") {
     
    qpcr_data <- qpcr_data[ , -1]
  
    }

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
colnames(qpcr_data)[match(ref1, colnames(qpcr_data))] <- "ref1"

if (ref2 != "") {

  colnames(qpcr_data)[match(ref2, colnames(qpcr_data))] <- "ref2"

}

colnames(qpcr_data)[match(goi, colnames(qpcr_data))] <- "goi"

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
control_ct <- qpcr_data[match(group1, qpcr_data$Sample), ]
qpcr_data <- dplyr::mutate(qpcr_data, avg_control_dCt = mean(control_ct$goi))

#ddCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, ddCt = dCt - avg_control_dCt)

#2^-ddCt calculation
qpcr_data <- dplyr::mutate(qpcr_data, expression = 2^-ddCt)

#mean the expression in biological replicates
control_exp <- qpcr_data[match(group1, qpcr_data$Sample), ]
target_exp1 <- qpcr_data[match(group2, qpcr_data$Sample), ]

if (group3 != "") {
  
 target_exp2 <- qpcr_data[match(group3, qpcr_data$Sample), ] 
  
}

if (group4 != "") {
  
 target_exp3 <- qpcr_data[match(group4, qpcr_data$Sample), ] 
  
}

if (group5 != "") {
  
 target_exp4 <- qpcr_data[match(group5, qpcr_data$Sample), ] 
  
}

control_exp <- dplyr::mutate(control_exp, control_avg_exp = mean(expression))
target_exp1 <- dplyr::mutate(target_exp1, target_avg_exp1 = mean(expression))

if (group3 != "") {

  target_exp2 <- dplyr::mutate(target_exp2, target_avg_exp2 = mean(expression))

}

if (group4 != "") {

  target_exp3 <- dplyr::mutate(target_exp3, target_avg_exp3 = mean(expression))

}

if (group5 != "") {

  target_exp4 <- dplyr::mutate(target_exp4, target_avg_exp4 = mean(expression))

}

qpcr_data <- dplyr::left_join(qpcr_data, control_exp)
qpcr_data <- dplyr::left_join(qpcr_data, target_exp1)

if (group3 != "") {

  qpcr_data <- dplyr::left_join(qpcr_data, target_exp2)

}

if (group4 != "") {

  qpcr_data <- dplyr::left_join(qpcr_data, target_exp3)

}

if (group5 != "") {

  qpcr_data <- dplyr::left_join(qpcr_data, target_exp4)

}

qpcr_data <- dplyr::mutate(qpcr_data, avg_exp = dplyr::coalesce(target_avg_exp1, control_avg_exp))

if (group3 != "") {

  qpcr_data <- dplyr::mutate(qpcr_data, avg_exp = dplyr::coalesce(target_avg_exp1,
								  target_avg_exp2,
								  control_avg_exp))

}

if (group4 != "") {

  qpcr_data <- dplyr::mutate(qpcr_data, avg_exp = dplyr::coalesce(target_avg_exp1,
								  target_avg_exp2,
								  target_avg_exp3,
								  control_avg_exp))

}

if (group5 != "") {

  qpcr_data <- dplyr::mutate(qpcr_data, avg_exp = dplyr::coalesce(target_avg_exp1,
								  target_avg_exp2,
								  target_avg_exp3,
								  target_avg_exp4,
								  control_avg_exp))

}

#converting raw expression to percentage for better visualization
qpcr_data <- dplyr::arrange(qpcr_data, match(qpcr_data$Sample, c(group1, group2, group3, group4, group5)))
qpcr_data <- dplyr::mutate(qpcr_data, percent_exp = qpcr_data$avg_exp/qpcr_data$control_avg_exp[1]*100)

if (test == TRUE) {

  #duplicating rows for accurate p-value calculation
  idx <- rep(1:nrow(target_exp1), tech_rep)
  target_exp1 <- target_exp1[idx, ]

  idx2 <- rep(1:nrow(control_exp), tech_rep)
  control_exp <- control_exp[idx2, ]

  if (stat == "t-test") {
    #here add pairwise t-test with an if statement for multiple comparisons

    #calculating p-value
    stats <- t.test(target_exp1$expression, control_exp$expression, alternative = "two.sided")

    #converting pvalues to *
    if (stats$p.value < 0.001) {

      pvalue <- "***"

      } else if (stats$p.value < 0.01) {

        pvalue <- "**"

        } else if (stats$p.value < 0.05) {

          pvalue <- "*"

          } else {

             pvalue <- "ns"

             }

  #this is the supported df layout for ggpubr::stat_pvalue_manuel()
  if (group3 == "" & group4 == "" & group5 == "") {

    qpcr_data2 <- tibble::tribble(~group1, ~group2, ~pvalue,
                                  group1, group2, pvalue)

    } else if (group4 == "" & group5 == "") {

      qpcr_data2 <- tibble::tribble(~group1, ~group2, ~group3, ~pvalue,
                                    group1, group2, group3, pvalue)

      } else if (group5 == "") {

        qpcr_data2 <- tibble::tribble(~group1, ~group2, ~group3, ~group4, ~pvalue,
                                      group1, group2, group3, group4, pvalue)

        }

  }

  #anova test here
  if (stat == "anova") {

  #first data wrangling for anova function


  #anova test
  stats <- aov()
  stats_anova <- TukeyHSD(stats)

  }

}

#calculating sd
target_sd1 <- sd(target_exp1$expression)
control_sd <- sd(control_exp$expression)

if (group3 != "") {

	target_sd2 <- sd(target_exp2$expression)
}

if (group4 != "") {

	target_sd3 <- sd(target_exp3$expression)
}

if (group5 != "") {

	target_sd4 <- sd(target_exp4$expression)

}

#distinct only 1 column
qpcr_data <- qpcr_data[!duplicated(qpcr_data$percent_exp), ]

#converting raw sd into percentage
qpcr_data <- tibble::column_to_rownames(qpcr_data, var = "Sample")
target_sd1 <- target_sd1/qpcr_data[group2, "avg_exp"]*100
control_sd <- control_sd/qpcr_data[group1, "avg_exp"]*100

if (ref2 != "") {

  percent_exp1 <- qpcr_data[group1, "percent_exp"]
  percent_exp2 <- qpcr_data[group2, "percent_exp"]

  if (group3 != "") {
    percent_exp3 <- qpcr_data[group3, "percent_exp"]
    }

  if (group4 != "") {
    percent_exp4 <- qpcr_data[group4, "percent_exp"]
    }

  if (group5 != "") {
    percent_exp5 <- qpcr_data[group5, "percent_exp"]
    }

  } else {

  percent_exp1 <- qpcr_data[group1, "percent_exp"]
  percent_exp2 <- qpcr_data[group2, "percent_exp"]

  if (group3 != "") {
    percent_exp3 <- qpcr_data[group3, "percent_exp"]
    }

  if (group4 != "") {
    percent_exp4 <- qpcr_data[group4, "percent_exp"]
    }

  if (group5 != "") {
    percent_exp5 <- qpcr_data[group5, "percent_exp"]
    }

}

qpcr_data <- tibble::rownames_to_column(qpcr_data, var = "Sample")

#removing na to only visualize sample of interest
qpcr_data <- qpcr_data[!is.na(qpcr_data$percent_exp), ]

#visualization
if (group3 == "" && group4 == "" && group5 == "") {

  if (test == TRUE) {

    final_plot <- ggpubr::ggbarplot(qpcr_data,
       	                            x = "Sample",
                                    y = "percent_exp",
			            fill = "808080",
				    xlab = "Sample",
				    ylab = "Relative mRNA level",
				    title = goi,
				    size = 0.5,
				    palette = "npg",
				    lab.size = 5,
				    lab.vjust = 0.5,
				    lab.hjust = 1.2,
				    sort.by.groups = FALSE,
				    ggtheme = ggpubr::theme_pubr(base_size = 14)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group2,
	        			 ymin = percent_exp2 - target_sd1,
	                		 ymax = percent_exp2 + target_sd1,
		                	 width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group1,
	    				 ymin = percent_exp1 - control_sd,
     					 ymax = percent_exp1 + control_sd,
     					 width = 0.1)) +
                  ggpubr::stat_pvalue_manual(qpcr_data2, label = "pvalue",
					     y.position = percent_exp2 + control_sd + 20,
					     bracket.size = 0.6, label.size = 6)

  } else {

    final_plot <- ggpubr::ggbarplot(qpcr_data,
				    x = "Sample",
				    y = "percent_exp",
				    fill = "808080",
				    xlab = "Sample",
				    ylab = "Relative mRNA level",
				    size = 0.5,
				    title = goi,
				    palette = "npg",
				    lab.size = 5,
				    lab.vjust = 0.5,
				    lab.hjust = 1.2,
				    sort.by.groups = FALSE,
				    ggtheme = ggpubr::theme_pubr(base_size = 14)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group2,
	        			 ymin = percent_exp2 - target_sd1,
	                		 ymax = percent_exp2 + target_sd1,
		                	 width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group1,
	    				 ymin = percent_exp1 - control_sd,
     					 ymax = percent_exp1 + control_sd,
     					 width = 0.1))

  }

} else if (group4 == "" && group5 == "") {

	if (test == TRUE) {

	} else {

    final_plot <- ggpubr::ggbarplot(qpcr_data,
				    x = "Sample",
				    y = "percent_exp",
				    fill = "808080",
				    xlab = "Sample",
				    ylab = "Relative mRNA level",
				    size = 0.5,
				    title = goi,
				    palette = "npg",
				    lab.size = 5,
				    lab.vjust = 0.5,
				    lab.hjust = 1.2,
				    sort.by.groups = FALSE,
				    ggtheme = ggpubr::theme_pubr(base_size = 14)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group2,
	        			              ymin = percent_exp2 - target_sd1,
	                		              ymax = percent_exp2 + target_sd1,
		                	              width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group3,
						      ymin = percent_exp3 - target_sd2,
						      ymax = percent_exp3 + target_sd2,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group1,
	    				              ymin = percent_exp1 - control_sd,
     					              ymax = percent_exp1 + control_sd,
     					              width = 0.1))

	}

} else if (group5 == "") {

	if (test == TRUE) {

	} else {

    final_plot <- ggpubr::ggbarplot(qpcr_data,
				    x = "Sample",
				    y = "percent_exp",
				    fill = "808080",
				    xlab = "Sample",
				    ylab = "Relative mRNA level",
				    title = goi,
				    size = 0.5,
				    palette = "npg",
				    lab.size = 5,
				    lab.vjust = 0.5,
				    lab.hjust = 1.2,
				    sort.by.groups = FALSE,
				    ggtheme = ggpubr::theme_pubr(base_size = 14)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group2,
	        			 ymin = percent_exp2 - target_sd1,
	                		 ymax = percent_exp2 + target_sd1,
		                	 width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group3,
						      ymin = percent_exp3 - target_sd2,
						      ymax = percent_exp3 + target_sd2,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group4,
						      ymin = percent_exp4 - target_sd3,
						      ymax = percent_exp4 + target_sd3,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group1,
	    				 ymin = percent_exp1 - control_sd,
     					 ymax = percent_exp1 + control_sd,
     					 width = 0.1))

	}

} else {

  if (test == TRUE) {

  } else {

    final_plot <- ggpubr::ggbarplot(qpcr_data,
				    x = "Sample",
				    y = "percent_exp",
				    fill = "808080",
				    xlab = "Sample",
				    ylab = "Relative mRNA level",
				    size = 0.5,
				    title = goi,
				    palette = "npg",
				    lab.size = 5,
				    lab.vjust = 0.5,
				    lab.hjust = 1.2,
				    sort.by.groups = FALSE,
				    ggtheme = ggpubr::theme_pubr(base_size = 14)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group2,
	        			 ymin = percent_exp2 - target_sd1,
	                		 ymax = percent_exp2 + target_sd1,
		                	 width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group3,
						      ymin = percent_exp3 - target_sd2,
						      ymax = percent_exp3 + target_sd2,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group3,
						      ymin = percent_exp4 - target_sd3,
						      ymax = percent_exp4 + target_sd3,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group3,
						      ymin = percent_exp5 - target_sd4,
						      ymax = percent_exp5 + target_sd4,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(ggplot2::aes(x = group1,
	    				 ymin = percent_exp1 - control_sd,
     					 ymax = percent_exp1 + control_sd,
     					 width = 0.1))

  }

}

plot(final_plot)
return(final_plot)

}
