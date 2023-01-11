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

allpops = read.pcadapt('Charr_Poly_All_Fixed_coords_maf05_geno95.bed', 
                           type = 'bed')

## Elbows at K=4
# pca_allpops = pcadapt(allpops, 
#                            K = 60, 
#                            method="mahalanobis", 
#                            min.maf = 0.01)
# plot(pca_allpops, 
#      option = 'screeplot')

pca_allpops = pcadapt(allpops, 
                      K = 10, 
                      method = 'mahalanobis', 
                      min.maf = 0.01)

plot(pca_allpops, 
     option = 'screeplot')

summary(pca_allpops)

pca_allpops$singular.values
sum(pca_allpops$singular.values)

(sqrt(0.51083013)/1.634673)*100
## 43.72%
(sqrt(0.23428876)/1.634673)*100
## 29.61%
(sqrt(0.19223258)/1.634673)*100
## 26.82%

pca_allpops_scores = as_tibble(pca_allpops$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_allpops_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_delim('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.ped')

meta_data = data %>% 
  dplyr::select(FID, 
                IndividualID) %>% 
  left_join(., 
            identifiers,
            by = '#FamilyID')

