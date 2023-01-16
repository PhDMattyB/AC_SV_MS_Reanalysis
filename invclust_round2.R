##############################
## Invclust across each chr
##
## Matt Brachmann (PhDMattyB)
##
## 2023-01-14
##
##############################

setwd('~/AC_SV_MS_Data/invclust')

library(patchwork)
library(invClust)
library(tidyverse)
library(snpStats)
BiocManager::install("inveRsion")
library(inveRsion)

theme_set(theme_bw())

## Not in function
'%!in%' <- function(x,y)!('%in%'(x,y))


# lg1 sv detect -----------------------------------------------------------


lg1_inv = read.plink(bed = 'Charr_Poly_All_Fixed_1.bed', 
                        bim = 'Charr_Poly_All_Fixed_1.bim', 
                        fam = 'Charr_Poly_All_Fixed_1.fam')

lg1_geno = lg1_inv$genotypes
lg1_map = lg1_inv$map 
## check
identical(lg1_map[,2],
          colnames(lg1_geno))

## 100 snps might be a lot to find any inversions

lg1_inv_detect = data.frame(chr = 1, 
                              LBP = 6793717, 
                              RBP = 84732430, 
                              reg = 'lg1_inversion')

lg1_inver1 = invClust(roi = lg1_inv_detect, 
                             wh = 1, 
                             geno = lg1_geno, 
                             annot = lg1_map, 
                             dim = 2)

plot(lg1_inver1)


