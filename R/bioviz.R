install.packages("devtools")
install.packages("roxygen2")
library(devtools)
library(roxygen2)
#biovizR is an easy to use publication/report ready
#biological data visualization package created by R.
#rather than relying on ggplot2 grammar.
#the package relies on functions such as
#qpcr_barplot(), rnaseq_maplot() etc etc.
#the purpose of this package is to provide easy to
#use R package aimed at wetlab biologists to explore and
#visualize their data
library(here)
install.packages("here")
devtools::create("biovizR")
here::dr_here()

