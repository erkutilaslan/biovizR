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
#' @param stat Default TRUE. Enable or disable statistical analysis.
#' @param test Default t-test. Select statistical testing method.
#' @param test_mode Default two.sided. Define alternative hypothesis for statistical test.
#' @param generate_table Default False. Set to True to create a table with the results everytime run the function.
#' @return A bar plot of qPCR results.
#' @import tidyverse
#' @export

#qpcr_data <- "~/N3_OE_results.csv"
#group1 <- "Control"
#group2 <- "NANOS3"
#group3 <- ""
#group4 <- ""
#group5 <- ""
#ref1 <- "GARS1"
#ref2 <- "DTD1"
#goi <- "FOXM1"
#stat <- TRUE
#test_mode <- "less"
#test <- "t-test"
#type <- "biorad"
#generate_table <- FALSE

#barplot_qpcr(qpcr_data, group1 = "Control", group2 = "NANOS3", ref1 = "GARS1", ref2 = "DTD1", goi = "FOXM1", stat = TRUE, test_mode = "less")

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
			 stat = TRUE,
                         test = "t-test",
			 test_mode = "two.sided",
			 generate_table = FALSE) {


####Functions####

#data import function
data_import <- function(import_data) {

	if (is.character(import_data) == TRUE) {
 
		if (grepl(".csv", as.character(import_data)) == TRUE) {
    
			import_data <- read.csv(import_data) 

		} else {

			import_data <- readxl::read_excel(import_data, sheet = 1)

		}

	}

	return(import_data)
}

#data wrangling function
data_process <- function(process_data) {

	if (type == "biorad") {

		if (colnames(process_data[1]) == "X" || colnames(process_data[1]) == "x") {
     
			process_data <- process_data[ , -1]
  
		}

		if (stat == TRUE) {

			colnames(process_data)[6] <- "Biological.Set.Name"
			process_data <- process_data[, c(3, 5, 6, 8)]
	
		} else {

			process_data <- process_data[, c(3, 5, 8)]

		}

		process_data$Sample[process_data$Sample == ""] <- NA
		process_data$Target[process_data$Target == "Target"] <- NA
		process_data <- na.omit(process_data)
		process_data$Cq.Mean <- as.numeric(as.character(process_data$Cq.Mean))

		if (ref2 == "") {

		goi_data <- dplyr::filter(process_data, Target == goi)
		ref1_data <- dplyr::filter(process_data, Target == ref1)
		process_data <- list(goi_data, ref1_data)
		
		} else {

		goi_data <- dplyr::filter(process_data, Target == goi)
		ref1_data <- dplyr::filter(process_data, Target == ref1)
		ref2_data <- dplyr::filter(process_data, Target == ref2)
		process_data <- list(goi_data, ref1_data, ref2_data)

		}
	
	}

	return(process_data)
}


ddCt_calc <- function (ddCt_data) {

	goi_data <- as.data.frame(ddCt_data[1])
	goi_data <- dplyr::arrange(goi_data, desc(Biological.Set.Name))
	goi_data <- dplyr::arrange(goi_data, desc(Sample))
	ref1_data <- as.data.frame(ddCt_data[2])
	ref1_data <- dplyr::arrange(ref1_data, desc(Biological.Set.Name))
	ref1_data <- dplyr::arrange(ref1_data, desc(Sample))

	#avg of refs for multiple refs
	if (ref2 != "") {

		ref2_data <- as.data.frame(ddCt_data[3])
		ref2_data <- dplyr::arrange(ref2_data, desc(Biological.Set.Name))
		ref2_data <- dplyr::arrange(ref2_data, desc(Sample))
		ddCt_data <- dplyr::mutate(goi_data, avg_ref = (ref1_data$Cq.Mean + ref2_data$Cq.Mean) / 2)
		#dCt calculation for avg of refs
		ddCt_data <- dplyr::mutate(ddCt_data, dCt = (Cq.Mean - avg_ref))

	} else {

		#dCt calculation for only 1 ref
		ddCt_data <- dplyr::mutate(goi_data, dCt = (goi_data$Cq.Mean - ref1_data$Cq.Mean))

	}

	#avg.dCt of target gene value in control biological replicates
	control_ct <- dplyr::filter(ddCt_data, Sample == group1)
	ddCt_data <- dplyr::mutate(ddCt_data, avg_control_dCt = mean(control_ct$Cq.Mean))
	ddCt_data <- dplyr::distinct(ddCt_data)

	#ddCt calculation
	ddCt_data <- dplyr::mutate(ddCt_data, ddCt = dCt - avg_control_dCt)

	#2^-ddCt calculation
	ddCt_data <- dplyr::mutate(ddCt_data, expression = 2^-ddCt)

	#mean the expression in biological replicates
	control_exp <- dplyr::filter(ddCt_data, Sample == group1)
	target_exp1 <- dplyr::filter(ddCt_data, Sample == group2)

	if (group3 != "") {
  
		target_exp2 <- dplyr::filter(ddCt_data, Sample == group3)
  
	}

	if (group4 != "") {
  
		target_exp3 <- dplyr::filter(ddCt_data, Sample == group4)
  
	}

	if (group5 != "") {
  
		target_exp4 <- dplyr::filter(ddCt_data, Sample == group5)
  
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

	ddCt_data <- dplyr::left_join(ddCt_data, control_exp)
	ddCt_data <- dplyr::left_join(ddCt_data, target_exp1)

	if (group3 != "") {

		ddCt_data <- dplyr::left_join(ddCt_data, target_exp2)

	}

	if (group4 != "") {

		ddCt_data <- dplyr::left_join(ddCt_data, target_exp3)

	}

	if (group5 != "") {

		ddCt_data <- dplyr::left_join(ddCt_data, target_exp4)

	}

	ddCt_data <- dplyr::mutate(ddCt_data, avg_exp = dplyr::coalesce(target_avg_exp1, control_avg_exp))

	if (group3 != "") {

		ddCt_data <- dplyr::mutate(ddCt_data, avg_exp = dplyr::coalesce(target_avg_exp1,
										target_avg_exp2,
										control_avg_exp))
	
	}

	if (group4 != "") {

		ddCt_data <- dplyr::mutate(ddCt_data, avg_exp = dplyr::coalesce(target_avg_exp1,
										target_avg_exp2,
										target_avg_exp3,
										control_avg_exp))

	}

	if (group5 != "") {

		ddCt_data <- dplyr::mutate(ddCt_data, avg_exp = dplyr::coalesce(target_avg_exp1,
									  	target_avg_exp2,
									  	target_avg_exp3,
									  	target_avg_exp4,
									  	control_avg_exp))

	}

	return(ddCt_data)
}

#raw expression to percentage function
raw_to_percent <- function(raw_to_percent_data) {


	if (group3 == "" & group4 == "" & group4 == "") {

		raw_to_percent_data <- dplyr::arrange(raw_to_percent_data,
						      match(raw_to_percent_data$Sample, c(group1, group2)))

	}

	if (group4 == "" & group5 == "") {

		raw_to_percent_data <- dplyr::arrange(raw_to_percent_data,
						      match(raw_to_percent_data$Sample, c(group1, group2, group3)))

	}

	if (group5 == "") {

		raw_to_percent_data <- dplyr::arrange(raw_to_percent_data,
						      match(raw_to_percent_data$Sample, c(group1, group2, group3, group4)))

	} else {

		raw_to_percent_data <- dplyr::arrange(raw_to_percent_data,
						      match(raw_to_percent_data$Sample, c(group1, group2, group3, group4, group5)))

	}
	percent_data_1 <- dplyr::filter(raw_to_percent_data, Biological.Set.Name == 1)
	percent_data_2 <- dplyr::filter(raw_to_percent_data, Biological.Set.Name == 2)
	percent_data_3 <- dplyr::filter(raw_to_percent_data, Biological.Set.Name == 3)

	

	percent_data_1 <- dplyr::mutate(percent_data_1,
					     percent_exp = percent_data_1$expression/percent_data_1$expression[1]*100)
	percent_data_2 <- dplyr::mutate(percent_data_2,
					     percent_exp = percent_data_2$expression/percent_data_2$expression[1]*100)
	percent_data_3 <- dplyr::mutate(percent_data_3,
					     percent_exp = percent_data_3$expression/percent_data_3$expression[1]*100)

	raw_to_percent_data <- dplyr::bind_rows(percent_data_1, percent_data_2, percent_data_3)
	raw_to_percent_data <- tidyr::drop_na(raw_to_percent_data, percent_exp)

	return(raw_to_percent_data)
}

#pvalue to star conversion function
pvalue_star <- function(dat) {

	if (dat < 0.001) {

		dat <- "***"

	} else if (dat < 0.01) {

		dat <- "**"

	} else if (dat < 0.05) {

		dat <- "*"

	} else if (dat >= 0.05){

		dat <- "ns"

	} else {

		dat <- NA

	}

}

#statistics function
qpcr_stat <- function(qpcr_stat_data) {

if (test == "t-test" & group3 == "" & group4 == "" & group5 == "") {

	control_exp <- dplyr::filter(qpcr_stat_data, Sample == group1)
	target_exp1 <- dplyr::filter(qpcr_stat_data, Sample == group2)

	#calculating p-value
	stats <- t.test(target_exp1$percent_exp, control_exp$percent_exp, alternative = test_mode, var.equal = TRUE)

	#converting pvalues to *
	pvalue <- pvalue_star(stats$p.value)

	#this is the supported df layout for ggpubr::stat_pvalue_manual()
	stat1 <- tibble::tribble(~group1, ~group2, ~pvalue,
						group1, group2, pvalue)
	stats <- stat1

	}

if (test == "t-test" & group3 != "" & group4 == "" & group5 == "") {

	target_exp2 <- dplyr::filter(qpcr_stat_data, Sample == group3)

	qpcr_data <- dplyr::full_join(qpcr_data, control_exp)
	qpcr_data <- dplyr::full_join(qpcr_data, target_exp1)
	qpcr_data <- dplyr::full_join(qpcr_data, target_exp2)
	stats <- pairwise.t.test(qpcr_data$expression, qpcr_data$Sample, p.adjust.method = "fdr")
	pvalues <- stats$p.value
	pvalues <- pvalues[,group1]

	pvalues <- lapply(pvalues, pvalue_star)
	pvalue1 <- as.character(pvalues[group2])
	pvalue2 <- as.character(pvalues[group3])

	stat1 <- tibble::tribble(~group1, ~group2, ~pvalue,
				 group1, group2, pvalue1)
	stat2 <- tibble::tribble(~group1, ~group2, ~pvalue,
				 group1, group3, pvalue2)

}

if (test == "t-test" & group3 != "" & group4 != "" & group5 == "") {

	target_exp2 <- dplyr::filter(qpcr_stat_data, Sample == group3)
	target_exp3 <- dplyr::filter(qpcr_stat_data, Sample == group4)

	qpcr_data <- dplyr::full_join(qpcr_data, control_exp)
	qpcr_data <- dplyr::full_join(qpcr_data, target_exp1)
	qpcr_data <- dplyr::full_join(qpcr_data, target_exp2)
	qpcr_data <- dplyr::full_join(qpcr_data, target_exp3)
	stats <- pairwise.t.test(qpcr_data$expression, qpcr_data$Sample, p.adjust.method = "fdr")
	pvalues <- stats$p.value
	pvalues <- pvalues[,group1]

	pvalues <- lapply(pvalues, pvalue_star)
	pvalue1 <- as.character(pvalues[group2])
	pvalue2 <- as.character(pvalues[group3])
	pvalue3 <- as.character(pvalues[group4])

	stat1 <- tibble::tribble(~group1, ~group2, ~pvalue,
				 group1, group2, pvalue1)
	stat2 <- tibble::tribble(~group1, ~group2, ~pvalue,
				 group1, group3, pvalue2)
	stat3 <- tibble::tribble(~group1, ~group2, ~pvalue,
				 group1, group4, pvalue3)

}

if (test == "t-test" & group3 != "" & group4 != "" & group5 != "") {

	target_exp2 <- dplyr::filter(qpcr_stat_data, Sample == group3)
	target_exp3 <- dplyr::filter(qpcr_stat_data, Sample == group4)
	target_exp4 <- dplyr::filter(qpcr_stat_data, Sample == group5)

	qpcr_data <- dplyr::full_join(qpcr_data, control_exp)
	qpcr_data <- dplyr::full_join(qpcr_data, target_exp1)
	qpcr_data <- dplyr::full_join(qpcr_data, target_exp2)
	qpcr_data <- dplyr::full_join(qpcr_data, target_exp3)
	qpcr_data <- dplyr::full_join(qpcr_data, target_exp4)
	stats <- pairwise.t.test(qpcr_data$expression, qpcr_data$Sample, p.adjust.method = "fdr")
	pvalues <- stats$p.value
	pvalues <- pvalues[,group1]

	pvalues <- lapply(pvalues, pvalue_star)
	pvalue1 <- as.character(pvalues[group2])
	pvalue2 <- as.character(pvalues[group3])
	pvalue3 <- as.character(pvalues[group4])
	pvalue4 <- as.character(pvalues[group5])

	stat1 <- tibble::tribble(~group1, ~group2, ~pvalue,
				group1, group2, pvalue1)

	stat2 <- tibble::tribble(~group1, ~group2, ~pvalue,
				group1, group3, pvalue2)

	stat3 <- tibble::tribble(~group1, ~group2, ~pvalue,
				group1, group4, pvalue3)

	stat4 <- tibble::tribble(~group1, ~group2, ~pvalue,
				group1, group5, pvalue4)

}

#anova test here
if (test == "anova") {

	#first data wrangling for anova function


	#anova test
	stats <- aov()
	stats_anova <- TukeyHSD(stats)

	}

return(stats)
}

#sd calculation function. this calculates SE not SD!!
calc_sd <- function(calc_sd_data) {

	control_exp <- dplyr::filter(calc_sd_data, Sample == group1)
	target_exp1 <- dplyr::filter(calc_sd_data, Sample == group2)

	control_sd <- sd(control_exp$expression)/sqrt(length(control_exp$expression))
	target_sd1 <- sd(target_exp1$expression)/sqrt(length(target_exp1$expression))
	control_exp <- dplyr::mutate(control_exp, avg_percent_exp = mean(percent_exp))
	target_exp1 <- dplyr::mutate(target_exp1, avg_percent_exp = mean(percent_exp))
	calc_sd_data <- dplyr::bind_rows(control_exp, target_exp1)


	if (group3 != "") {

		target_exp2 <- dplyr::filter(calc_sd_data, Sample == group3)
		target_sd2 <- sd(target_exp2$expression)

	}

	if (group4 != "") {

		target_exp2 <- dplyr::filter(calc_sd_data, Sample == group4)
		target_sd3 <- sd(target_exp3$expression)

	}

	if (group5 != "") {

		target_exp4 <- dplyr::filter(calc_sd_data, Sample == group5)
		target_sd4 <- sd(target_exp4$expression)

	}

	#distinct only 1 column
	
	calc_sd_data <- calc_sd_data[!duplicated(calc_sd_data$avg_percent_exp), ]

	#converting raw sd into percentage
	rownames(calc_sd_data) <- NULL
	calc_sd_data <- tibble::column_to_rownames(calc_sd_data, var = "Sample")
	percent_sd1 <- control_sd*calc_sd_data[group1, "percent_exp"]/calc_sd_data[group1, "avg_exp"]
	percent_sd2 <- target_sd1*calc_sd_data[group2, "percent_exp"]/calc_sd_data[group2, "avg_exp"]

	if (group3 == "" & group4 == "" & group5 == "") {

		percent_sd <- c(percent_sd1, percent_sd2)

	}

	if (group3 != "") {
	
		percent_sd3 <- target_sd2/calc_sd_data[group3, "percent_exp"]/calc_sd_data[group3, "avg_exp"]
		percent_sd <- c(percent_sd1, percent_sd2, percent_sd3)

	}

	if (group4 != "") {

		percent_sd4 <- target_sd3/calc_sd_data[group4, "percent_exp"]/calc_sd_data[group4, "avg_exp"]
		percent_sd <- c(percent_sd1, percent_sd2, percent_sd3, percent_sd4)

	}

	if (group5 != "") {

		percent_sd5 <- target_sd4/calc_sd_data[group5, "percent_exp"]/calc_sd_data[group5, "avg_exp"]
		percent_sd <- c(percent_sd1, percent_sd2, percent_sd3, percent_sd4, percent_sd5)

	}
	

	return(percent_sd)
}

#get exp data as a vector from qpcr_data. needed for visualization
get_exp <- function(get_exp_data) {

	control_exp <- dplyr::filter(get_exp_data, Sample == group1)
	target_exp1 <- dplyr::filter(get_exp_data, Sample == group2)
	control_exp <- dplyr::mutate(control_exp, avg_percent_exp = mean(percent_exp))
	target_exp1 <- dplyr::mutate(target_exp1, avg_percent_exp = mean(percent_exp))
	get_exp_data <- dplyr::bind_rows(control_exp, target_exp1)
	#distinct only 1 column
	get_exp_data <- get_exp_data[!duplicated(get_exp_data$avg_percent_exp), ]
	rownames(get_exp_data) <- NULL
	get_exp_data <- tibble::column_to_rownames(get_exp_data, var = "Sample")

	percent_exp1 <- get_exp_data[group1, "avg_percent_exp"]
	percent_exp2 <- get_exp_data[group2, "avg_percent_exp"]
	
	if (group3 == "" & group4 == "" & group5 == "") {

		percent_exp <- c(percent_exp1, percent_exp2)

	}

	if (group3 != "") {

		percent_exp3 <- qpcr_data[group3, "percent_exp"]
		percent_exp <- c(percent_exp1, percent_exp2, percent_exp3)

	}

	if (group4 != "") {

   		percent_exp4 <- qpcr_data[group4, "percent_exp"]
		percent_exp <- c(percent_exp1, percent_exp2, percent_exp3, percent_exp4)

	}

	if (group5 != "") {

		percent_exp5 <- qpcr_data[group5, "percent_exp"]
		percent_exp <- c(percent_exp1, percent_exp2, percent_exp3, percent_exp4, percent_exp5)

	}

	get_exp_data <- tibble::rownames_to_column(get_exp_data, var = "Sample")

	return(percent_exp)
}

#import table function
export_table <- function(export_data, export_stat, export_sd) {

	export_data <- export_data[!duplicated(export_data$percent_exp), ]
	#generating results as a table
	samples <- as.character(export_data$Sample)
	expressions <- export_data$percent_exp
	sds <- export_sd
	stats <- c(export_stat$pvalue, export_stat$pvalue)
	targets <- c(goi, goi)
	results_table <- data.frame(samples, expressions, sds, stats, targets)
	write.table(results_table, quote = FALSE, sep = ",", file = paste0(goi, "_table.csv"), row.names = FALSE)

}

#plotting function
qpcr_plot <- function(qpcr_plot_data, qpcr_plot_exp, qpcr_plot_stat, qpcr_plot_sd) {

	control_exp <- dplyr::filter(qpcr_plot_data, Sample == group1)
	target_exp1 <- dplyr::filter(qpcr_plot_data, Sample == group2)
	control_exp <- dplyr::mutate(control_exp, avg_percent_exp = mean(percent_exp))
	target_exp1 <- dplyr::mutate(target_exp1, avg_percent_exp = mean(percent_exp))
	qpcr_plot_data <- dplyr::bind_rows(control_exp, target_exp1)
	#distinct only 1 column
	qpcr_plot_data <- qpcr_plot_data[!duplicated(qpcr_plot_data$avg_percent_exp), ]

if (group3 == "" && group4 == "" && group5 == "") {

	if (stat == TRUE) {

		final_plot <- ggpubr::ggbarplot(qpcr_plot_data,
       						x = "Sample",
						y = "avg_percent_exp",
						fill = "808080",
						xlab = "Sample",
						ylab = "Relative mRNA level",
						title = goi,
						size = 1,
						sort.by.groups = FALSE,
						ggtheme = ggpubr::theme_pubr(base_size = 18)) +
                ggplot2::theme(axis.line =  ggplot2::element_line(size = 1)) +
		ggplot2::geom_errorbar(size = 1, ggplot2::aes(x = group2,
						    ymin = qpcr_plot_exp[2] - qpcr_plot_sd[2],
						    ymax = qpcr_plot_exp[2] + qpcr_plot_sd[2],
						    width = 0.1)) +
		ggplot2::geom_errorbar(size = 1, ggplot2::aes(x = group1,
						    ymin = qpcr_plot_exp[1] - qpcr_plot_sd[1],
						    ymax = qpcr_plot_exp[1] + qpcr_plot_sd[1],
						    width = 0.1)) +
		ggpubr::stat_pvalue_manual(qpcr_plot_stat, label = "pvalue",
					   y.position = max(qpcr_plot_exp) + max(percent_sd) + 10,
					   bracket.size = 1, label.size = 8)

	} else {

		final_plot <- ggpubr::ggbarplot(qpcr_plot_data,
						x = "Sample",
						y = "avg_percent_exp",
						fill = "808080",
						xlab = "Sample",
						ylab = "Relative mRNA level",
						size = 1,
						title = goi,
						palette = "npg",
						lab.size = 5,
						lab.vjust = 0.5,
						lab.hjust = 1.2,
						sort.by.groups = FALSE,
						ggtheme = ggpubr::theme_pubr(base_size = 18)) +
                ggplot2::theme(axis.line =  ggplot2::element_line(size = 1)) +
		ggplot2::geom_errorbar(size = 1,
				       ggplot2::aes(x = group2,
	        				    ymin = qpcr_plot_exp[2] - qpcr_plot_sd[2],
						    ymax = qpcr_plot_exp[2] + qpcr_plot_sd[2],
						    width = 0.1)) +
		ggplot2::geom_errorbar(size = 1,
				       ggplot2::aes(x = group1,
		    				    ymin = qpcr_plot_exp[1] - qpcr_plot_sd[1],
     					 	    ymax = qpcr_plot_exp[1] + qpcr_plot_sd[1],
     					 	    width = 0.1))

	}

} else if (group4 == "" && group5 == "") {

	if (stat == TRUE) {

	final_plot <- ggpubr::ggbarplot(qpcr_plot_data,
					x = "Sample",
					y = "avg_percent_exp",
					fill = "808080",
					xlab = "Sample",
					ylab = "Relative mRNA level",
					size = 1,
					title = goi,
					palette = "npg",
					lab.size = 5,
					lab.vjust = 0.5,
					lab.hjust = 1.2,
					sort.by.groups = FALSE,
					ggtheme = ggpubr::theme_pubr(base_size = 18)) +
                ggplot2::theme(axis.line =  ggplot2::element_line(size = 1)) +
			ggplot2::geom_errorbar(size = 1,
					       ggplot2::aes(x = group2,
		        				    ymin = qpcr_plot_exp[2] - qpcr_plot_sd[2],
	                				    ymax = qpcr_plot_exp[2] + qpcr_plot_sd[2],
							    width = 0.1)) +
			ggplot2::geom_errorbar(size = 1,
					       ggplot2::aes(x = group3,
						    	    ymin = percent_exp3 - percent_sd3,
							    ymax = percent_exp3 + percent_sd3,
							    width = 0.1)) +
			ggplot2::geom_errorbar(size = 1,
					       ggplot2::aes(x = group1,
	    						    ymin = qpcr_plot_exp[1] - qpcr_plot_sd[1],
     						       	    ymax = qpcr_plot_exp[1] + qpcr_plot_sd[1],
     					                    width = 0.1)) +
			ggpubr::stat_pvalue_manual(qpcr_plot_stat, label = "pvalue",
						   y.position = max(qpcr_plot_exp) + max(percent_sd) + 10,
						   bracket.size = 1, label.size = 8)
			ggpubr::stat_pvalue_manual(stat2, label = "pvalue",
						   y.position = max(qpcr_plot_exp) + max(percent_sd) + 20,
						   bracket.size = 1, label.size = 8) 

	} else {

	final_plot <- ggpubr::ggbarplot(qpcr_plot_data,
					x = "Sample",
					y = "avg_percent_exp",
					fill = "808080",
					xlab = "Sample",
					ylab = "Relative mRNA level",
					size = 1,
					title = goi,
					palette = "npg",
					lab.size = 5,
					lab.vjust = 0.5,
					lab.hjust = 1.2,
					sort.by.groups = FALSE,
					ggtheme = ggpubr::theme_pubr(base_size = 18)) +
                ggplot2::theme(axis.line =  ggplot2::element_line(size = 1)) +
			ggplot2::geom_errorbar(size = 1,
					       ggplot2::aes(x = group2,
		        				    ymin = qpcr_plot_exp[2] - qpcr_plot_sd[2],
	                				    ymax = qpcr_plot_exp[2] + qpcr_plot_sd[2],
							    width = 0.1)) +
			ggplot2::geom_errorbar(size = 1,
					       ggplot2::aes(x = group3,
						    	    ymin = percent_exp3 - percent_sd3,
							    ymax = percent_exp3 + percent_sd3,
							    width = 0.1)) +
			ggplot2::geom_errorbar(size = 1,
					       ggplot2::aes(x = group1,
	    						    ymin = qpcr_plot_exp[1] - qpcr_plot_sd[1],
     						       	    ymax = qpcr_plot_exp[1] + qpcr_plot_sd[1],
     					                    width = 0.1))

	}

} else if (group5 == "") {

	if (stat == TRUE) {

		final_plot <- ggpubr::ggbarplot(qpcr_plot_data,
						x = "Sample",
						y = "avg_percent_exp",
						fill = "808080",
						xlab = "Sample",
						ylab = "Relative mRNA level",
						title = goi,
						size = 1,
						palette = "npg",
						lab.size = 5,
						lab.vjust = 0.5,
						lab.hjust = 1.2,
						sort.by.groups = FALSE,
						ggtheme = ggpubr::theme_pubr(base_size = 18)) +
                ggplot2::theme(axis.line =  ggplot2::element_line(size = 1)) +
				ggplot2::geom_errorbar(size = 1,
						       ggplot2::aes(x = group2,
		        			   		    ymin = qpcr_plot_exp[2] - qpcr_plot_sd[2],
		                		 		    ymax = qpcr_plot_exp[2] + qpcr_plot_sd[2],
					                	    width = 0.1)) +
				ggplot2::geom_errorbar(size = 1,
						       ggplot2::aes(x = group3,
							      	    ymin = percent_exp3 - percent_sd3,
							      	    ymax = percent_exp3 + percent_sd3,
							      	    width = 0.1)) +
				ggplot2::geom_errorbar(size = 1,
						       ggplot2::aes(x = group4,
							      	    ymin = percent_exp4 - percent_sd4,
							      	    ymax = percent_exp4 + percent_sd4,
							      	    width = 0.1)) +
				ggplot2::geom_errorbar(size = 1,
						       ggplot2::aes(x = group1,
		    				 		    ymin = qpcr_plot_exp[1] - qpcr_plot_sd[1],
     						 		    ymax = qpcr_plot_exp[1] + qpcr_plot_sd[1],
     						 		    width = 0.1)) +
				ggpubr::stat_pvalue_manual(qpcr_plot_stat, label = "pvalue",
						           y.position = max(qpcr_plot_exp) + max(percent_sd) + 10,
					   bracket.size = 1, label.size = 8) +
				ggpubr::stat_pvalue_manual(stat2, label = "pvalue",
						           y.position = max(qpcr_plot_exp) + max(percent_sd) + 20,
					   bracket.size = 1, label.size = 8) +
				ggpubr::stat_pvalue_manual(stat3, label = "pvalue",
						           y.position = max(qpcr_plot_exp) + max(percent_sd) + 30,
							   bracket.size = 1, label.size = 8)

	} else {

		final_plot <- ggpubr::ggbarplot(qpcr_plot_data,
						x = "Sample",
						y = "avg_percent_exp",
						fill = "808080",
						xlab = "Sample",
						ylab = "Relative mRNA level",
						title = goi,
						size = 1,
						palette = "npg",
						lab.size = 5,
						lab.vjust = 0.5,
						lab.hjust = 1.2,
						sort.by.groups = FALSE,
						ggtheme = ggpubr::theme_pubr(base_size = 18)) +
                ggplot2::theme(axis.line =  ggplot2::element_line(size = 1)) +
				ggplot2::geom_errorbar(size = 1,
						       ggplot2::aes(x = group2,
		        			   		    ymin = qpcr_plot_exp[2] - qpcr_plot_sd[2],
		                		 		    ymax = qpcr_plot_exp[2] + qpcr_plot_sd[2],
					                	    width = 0.1)) +
				ggplot2::geom_errorbar(size = 1,
						       ggplot2::aes(x = group3,
							      	    ymin = percent_exp3 - percent_sd3,
							      	    ymax = percent_exp3 + percent_sd3,
							      	    width = 0.1)) +
				ggplot2::geom_errorbar(size = 1,
						       ggplot2::aes(x = group4,
							      	    ymin = percent_exp4 - percent_sd4,
							      	    ymax = percent_exp4 + percent_sd4,
							      	    width = 0.1)) +
				ggplot2::geom_errorbar(size = 1,
						       ggplot2::aes(x = group1,
		    				 		    ymin = qpcr_plot_exp[1] - qpcr_plot_sd[1],
     						 		    ymax = qpcr_plot_exp[1] + qpcr_plot_sd[1],
     						 		    width = 0.1))

	}

} else {
			
	if (stat == TRUE) {

		final_plot <- ggpubr::ggbarplot(qpcr_plot_data,
						x = "Sample",
				 		y = "avg_percent_exp",
						fill = "808080",
				    		xlab = "Sample",
				    		ylab = "Relative mRNA level",
				    		size = 1,
				    		title = goi,
				    		palette = "npg",
				    		lab.size = 5,
				    		lab.vjust = 0.5,
				    		lab.hjust = 1.2,
				    		sort.by.groups = FALSE,
				    		ggtheme = ggpubr::theme_pubr(base_size = 18)) +
                ggplot2::theme(axis.line =  ggplot2::element_line(size = 1)) +
                  ggplot2::geom_errorbar(size = 1,
					 ggplot2::aes(x = group2,
	        			              ymin = qpcr_plot_exp[2] - qpcr_plot_sd[2],
	                		 	      ymax = qpcr_plot_exp[2] + qpcr_plot_sd[2],
		                	 	      width = 0.1)) +
                  ggplot2::geom_errorbar(size = 1,
					 ggplot2::aes(x = group3,
						      ymin = percent_exp3 - percent_sd3,
						      ymax = percent_exp3 + percent_sd3,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(size = 1,
					 ggplot2::aes(x = group4,
						      ymin = percent_exp4 - percent_sd4,
						      ymax = percent_exp4 + percent_sd4,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(size = 1,
					 ggplot2::aes(x = group5,
						      ymin = percent_exp5 - percent_sd5,
						      ymax = percent_exp5 + percent_sd5,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(size = 1,
					 ggplot2::aes(x = group1,
	    				 	      ymin = qpcr_plot_exp[1] - qpcr_plot_sd[1],
     					 	      ymax = qpcr_plot_exp[1] + qpcr_plot_sd[1],
     					 	      width = 0.1)) +
		ggpubr::stat_pvalue_manual(qpcr_plot_stat, label = "pvalue",
				           y.position = max(qpcr_plot_exp) + max(percent_sd) + 10,
					   bracket.size = 1, label.size = 8) +
		ggpubr::stat_pvalue_manual(stat2, label = "pvalue",
				           y.position = max(qpcr_plot_exp) + max(percent_sd) + 20,
					   bracket.size = 1, label.size = 8) +
		ggpubr::stat_pvalue_manual(stat3, label = "pvalue",
				           y.position = max(qpcr_plot_exp) + max(percent_sd) + 30,
					   bracket.size = 1, label.size = 8) +
		ggpubr::stat_pvalue_manual(stat4, label = "pvalue",
				           y.position = max(qpcr_plot_exp) + max(percent_sd) + 40,
					   bracket.size = 1, label.size = 8)

	} else {

		final_plot <- ggpubr::ggbarplot(qpcr_plot_data,
						x = "Sample",
				 		y = "avg_percent_exp",
						fill = "808080",
				    		xlab = "Sample",
				    		ylab = "Relative mRNA level",
				    		size = 1,
				    		title = goi,
				    		palette = "npg",
				    		lab.size = 5,
				    		lab.vjust = 0.5,
				    		lab.hjust = 1.2,
				    		sort.by.groups = FALSE,
				    		ggtheme = ggpubr::theme_pubr(base_size = 18)) +
                ggplot2::theme(axis.line =  ggplot2::element_line(size = 1)) +
                  ggplot2::geom_errorbar(size = 1,
					 ggplot2::aes(x = group2,
	        			              ymin = qpcr_plot_exp[2] - qpcr_plot_sd[2],
	                		 	      ymax = qpcr_plot_exp[2] + qpcr_plot_sd[2],
						      width = 0.1)) +
                  gplot2::geom_errorbar(size = 1,
					ggplot2::aes(x = group3,
						      ymin = percent_exp3 - percent_sd3,
						      ymax = percent_exp3 + percent_sd3,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(size = 1,
					 ggplot2::aes(x = group4,
						      ymin = percent_exp4 - percent_sd4,
						      ymax = percent_exp4 + percent_sd4,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(size = 1,
					 ggplot2::aes(x = group5,
						      ymin = percent_exp5 - percent_sd5,
						      ymax = percent_exp5 + percent_sd5,
						      width = 0.1)) +
                  ggplot2::geom_errorbar(size = 1,
					 ggplot2::aes(x = group1,
	    				 	      ymin = qpcr_plot_exp[1] - qpcr_plot_sd[1],
     					 	      ymax = qpcr_plot_exp[1] + qpcr_plot_sd[1],
     					 	      width = 0.1))

	}

}

return(final_plot)
}


#test pipeline
qpcr_data <- data_import(qpcr_data)
qpcr_data <- data_process(qpcr_data)
qpcr_data <- ddCt_calc(qpcr_data)
qpcr_data <- raw_to_percent(qpcr_data)
stats <- qpcr_stat(qpcr_data)
percent_sd <- calc_sd(qpcr_data)
exp_data <- get_exp(qpcr_data)
final_plot_qpcr <- qpcr_plot(qpcr_data, exp_data, stats, percent_sd)

if (generate_table == TRUE) {
export_table(qpcr_data, stats, percent_sd)
}
plot(final_plot_qpcr)
return(final_plot_qpcr)

}
