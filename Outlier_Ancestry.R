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



# Fst setup outlier loci filter -------------------------------------------

read_csv('AC_bioclim_partial_RDA_outlier_data_25.06.2022.csv') %>% 
  dplyr::select(Chromosome, 
                SNP...1, 
                Genetic_pos, 
                Position) %>% 
  rename(SNP = SNP...1, 
         BP = Position) %>% 
  arrange(Chromosome) %>% 
  dplyr::select(SNP) %>% 
  write_tsv('bioclim_outlier_keep.txt')
  

read_csv('AC_sdm_partial_RDA_outlier_data_25.01.2023.csv') %>% 
  dplyr::select(Chromosome, 
                SNP...1, 
                Genetic_pos, 
                Position) %>% 
  rename(SNP = SNP...1, 
         BP = Position) %>% 
  arrange(Chromosome) %>%
  dplyr::select(SNP) %>% 
  write_tsv('sdm_outlier_keep.txt')



# Genomewide ancestry -----------------------------------------------------




# bioclim outlier ancestry ------------------------------------------------




# sdm outlier ancestry ----------------------------------------------------


