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
ATL_ACD  = read_tsv('ATL_ACD_Fst.fst') %>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

ARC_ACD  = read_tsv('ARC_ACD_Fst.fst') %>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

ARC_ATL  = read_tsv('ARC_ATL_Fst.fst') %>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

# sliding window analysis -------------------------------------------------


fst_100Kb = winScan(x = ARC_ATL, 
                          groups = 'CHR', 
                          position = 'POS',
                          values = 'FST_zero', 
                          win_size = 100000, 
                          win_step = 99999, 
                          funs = c('mean', 'sd'))

fst_100Kb = fst_100Kb %>%
  as_tibble() %>% 
  filter(FST_zero_n >= 3) %>% 
  write_tsv('ARC_ATL_Fst_100kb.txt')



# Fst violin plot ---------------------------------------------------------

ATL_ACD = read_tsv('ATL_ACD_Fst_100kb.txt')
ARC_ACD = read_tsv('ARC_ACD_Fst_100kb.txt')
ARC_ATL = read_tsv('ARC_ATL_Fst_100kb.txt')


