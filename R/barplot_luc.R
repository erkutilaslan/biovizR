#' This function visualizes qPCR results as bar plot.
#'
#' @param luc_data Input qPCR data to visualize.
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
#' @param stat Default TRUE. Enable or disable statistical analysis.
#' @param test Default t-test. Select statistical testing method.
#' @param generate_table Default False. Set to True to create a table with the results everytime run the function.
#' @return A bar plot of qPCR results.
#' @import tidyverse
#' @export

luc_data <- read.csv("C://Users/Erkut Ilaslan/Desktop/testluc_Dual-Luciferase 2 injectors_8-13-2021_12-07-26 PM - Copy.csv")

barplot_luc <- function(luc_data,
                        typle = "glomax",
                        group1 = "",
                        group2 = "",
                        group3 = "",
                        group4 = "",
                        group5 = "",
                        stat = TRUE,
                        test = "t-test",
                        generate_table = FALSE) {
  
#data import
  
  if (is.character(luc_data) == TRUE) {
    
    if (grepl(".csv", as.character(luc_data)) == TRUE) {
      
      luc_data <- read.csv(luc_data) 
      
    } else {
      
      luc_data <- readxl::read_excel(luc_data, sheet = 1)
      
    }
    
  }

#data wrangling
luc_data <- luc_data[33:41, 1:12]
luc_data <- luc_data[-1, -1:-2]
luc_data <- luc_data[luc_data == "X"] <- NA
luc_data <- luc_data[complete.cases(luc_data), ]
  
  
  
  
}
