##############################
##  Glacial lineage Fst
##
## Matt Brachmann (PhDMattyB)
##
## 14.04.2023
##
##############################

setwd('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/')

library(tidyverse)


# Make Keep files for plink -----------------------------------------------

ped_ids = read_table2('Charr_Poly_All_Fixed_coords_maf05_geno95.fam', 
                      col_names = F) %>%
  dplyr::select(X1,
                X2) %>% 
  rename(Population = X1, 
         Individual = X2)

latlong = read_csv('~/Bradbury_Postdoc/AC_SV_MS_Data/Admixture/SampleSites_Coords_1June2020.csv')


Data = inner_join(ped_ids, 
                  latlong) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin)

Data %>% 
  filter(Glacial_lin %in% c('ATL', 
                            'ACD')) %>% 
  dplyr::select(Population, 
                Individual) %>% 
  write_tsv('ATL_ACD_keep.txt', 
            col_names = F)

Data %>% 
  filter(Glacial_lin %in% c('ARC', 
                            'ACD')) %>% 
  dplyr::select(Population, 
                Individual) %>% 
  write_tsv('ARC_ACD_keep.txt', 
            col_names = F)

Data %>% 
  filter(Glacial_lin %in% c('ARC', 
                            'ATL')) %>% 
  dplyr::select(Population, 
                Individual) %>% 
  write_tsv('ARC_ATL_keep.txt', 
            col_names = F)




# Make Fst files for plink ------------------------------------------------

Data %>% 
  filter(Glacial_lin %in% c('ATL', 
                            'ACD')) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin) %>% 
  write_tsv('ATL_ACD_Fst_Grouping.txt', 
            col_names = F)

Data %>% 
  filter(Glacial_lin %in% c('ARC', 
                            'ACD')) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin) %>% 
  write_tsv('ARC_ACD_Fst_Grouping.txt', 
            col_names = F)

Data %>% 
  filter(Glacial_lin %in% c('ARC', 
                            'ATL')) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin) %>% 
  write_tsv('ARC_ATL_Fst_Grouping.txt', 
            col_names = F)
 


# Fst data ----------------------------------------------------------------




# sliding window analysis -------------------------------------------------



