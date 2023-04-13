##############################
## Outlier ancestry 
##
## Matt Brachmann (PhDMattyB)
##
## 13.04.2023
##
##############################

# Make --keep files for plink -------------------------------

# ped_test = read_table2('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.ped', 
#                        col_names = F)

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
  filter(Glacial_lin %in% c('Hybrid', 
                            'ATL')) %>% 
  dplyr::select(Population, 
                Individual) %>% 
  write_tsv('Lab_ATL_keep.txt', 
            col_names = F)

Data %>% 
  filter(Glacial_lin %in% c('Hybrid', 
                            'ARC')) %>% 
  dplyr::select(Population, 
                Individual) %>% 
  write_tsv('Lab_ARC_keep.txt', 
            col_names = F)

Data %>% 
  filter(Glacial_lin %in% c('Hybrid', 
                            'ACD')) %>% 
  dplyr::select(Population, 
                Individual) %>% 
  write_tsv('Lab_ACD_keep.txt', 
            col_names = F)


# Fst set up Plink ---------------------------------------------------------------

Data %>% 
  filter(Glacial_lin %in% c('Hybrid', 
                            'ATL')) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin) %>% 
  write_tsv('Lab_ATL_Fst_Grouping.txt', 
            col_names = F)

Data %>% 
  filter(Glacial_lin %in% c('Hybrid', 
                            'ARC')) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin) %>% 
  write_tsv('Lab_ARC_Fst_Grouping.txt', 
            col_names = F)

Data %>% 
  filter(Glacial_lin %in% c('Hybrid', 
                            'ACD')) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin) %>% 
  write_tsv('Lab_ACD_Fst_Grouping.txt', 
            col_names = F)
