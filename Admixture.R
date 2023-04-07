##############################
## Admixture results
##
## Matt Brachmann (PhDMattyB)
##
## 06.04.2023
##
##############################


setwd('~/Bradbury_Postdoc/AC_SV_MS_Data/Admixture/')


library(tidyverse)


qvalues = read_table('Charr_Poly_All_Fixed_coords_maf05_geno95.4.Q', 
                   col_names = c('Q1', 
                                 'Q2', 
                                 'Q3', 
                                 'Q4'))


