##############################
## Putative structural variant reanalysis
##
## Matt Brachmann (PhDMattyB)
##
## 2023-01-10
##
##############################

setwd('~/AC_SV_MS_Data/Pcadapt')

library(data.table)
library(tidyverse)
library(pcadapt)
library(devtools)
# install_github('isglobal-brge/invClust')
library(invClust)

theme_set(theme_bw())



# All populations pca -----------------------------------------------------

pca_allpops = read.pcadapt('Charr_Poly_All_Fixed_coords_maf05_geno95.bed', 
                           type = 'bed')

