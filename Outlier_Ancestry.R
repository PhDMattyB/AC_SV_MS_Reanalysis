##############################
## Outlier ancestry 
##
## Matt Brachmann (PhDMattyB)
##
## 13.04.2023
##
##############################

library(tidyverse)
library(windowscanr)

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



# Genomewide fst sliding window -----------------------------------------------------

lab_atl_fst = read_tsv('Lab_ATL_Fst.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

lab_arc_fst = read_tsv('Lab_ARC_Fst.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

lab_acd_fst = read_tsv('Lab_ACd_Fst.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))



fst_1kb = winScan(x = lab_atl_fst, 
                          groups = 'CHR', 
                          position = 'POS',
                          values = 'FST_zero', 
                          win_size = 1000, 
                          win_step = 999, 
                          funs = c('mean', 'sd'))

fst_1kb = fst_1kb %>%
  as_tibble() %>% 
  filter(FST_n >= 3) %>% 
  write_tsv('lab_atl_genomewide_fst_1kb.txt')



# Genomewide fst  -----------------------------------------------------

lab_atl_fst = read_tsv('Lab_ATL_Fst.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

lab_arc_fst = read_tsv('Lab_ARC_Fst.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

lab_acd_fst = read_tsv('Lab_ACd_Fst.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

genomewide_fst = ggplot() +
  geom_density(data = lab_atl_fst, 
               aes(x = FST_zero), 
               fill="#219ebc",
               alpha=0.8)+
  geom_density(data = lab_arc_fst, 
               aes(x = FST_zero), 
               fill = '#ff006e', 
               alpha = 0.8)+
  geom_density(data = lab_acd_fst, 
               aes(x = FST_zero), 
               fill = '#fb8500', 
               alpha = 0.8)+
  labs(x = 'Per locus Fst', 
       y = 'Density', 
       title = 'A)')+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())


# bioclim outlier ancestry ------------------------------------------------
lab_atl_bioclim = read_tsv('Lab_ATL_Fst_bioclim_outs.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

lab_arc_bioclim = read_tsv('Lab_ARC_Fst_bioclim_outs.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

lab_acd_bioclim = read_tsv('Lab_ACD_Fst_bioclim_outs.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

theme_set(theme_bw())

bioclim_outlier_fst = ggplot(data = lab_atl_bioclim, 
       aes(x = FST_zero)) +
  geom_density(fill="#219ebc",
               alpha=0.8)+
  geom_density(data = lab_arc_bioclim, 
              aes(x = FST_zero), 
              fill = '#ff006e', 
              alpha = 0.8)+
  geom_density(data = lab_acd_bioclim, 
               aes(x = FST_zero), 
               fill = '#fb8500', 
               alpha = 0.8)+
  labs(x = 'Per locus Fst', 
       y = 'Density', 
       title = 'B)')+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())



##

# sdm outlier ancestry ----------------------------------------------------

lab_atl_sdm = read_tsv('Lab_ATL_Fst_sdm_outs.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

lab_arc_sdm = read_tsv('Lab_ARC_Fst_sdm_outs.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

lab_acd_sdm = read_tsv('Lab_ACD_Fst_sdm_outs.fst') %>% 
  # na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))

theme_set(theme_bw())

sdm_outlier_fst = 
  ggplot() +
  geom_density(data = lab_arc_sdm, 
               aes(x = FST_zero), 
               fill = '#ff006e', 
               alpha = 0.8)+
  geom_density(data = lab_atl_sdm, 
               aes(x = FST_zero), 
               fill="#219ebc",
               alpha=0.8)+
  geom_density(data = lab_acd_sdm, 
               aes(x = FST_zero), 
               fill = '#fb8500', 
               alpha = 0.8)+
  labs(x = 'Per locus Fst', 
       y = 'Density', 
       title = 'C)')+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14))


# combine plots -----------------------------------------------------------

fst_plots = genomewide_fst/bioclim_outlier_fst/sdm_outlier_fst  


ggsave('Fst_patterns_genomewide_outliers.tiff', 
       plot = fst_plots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 15, 
       height = 15)

