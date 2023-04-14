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
label1 = rep('ATL-ACD', nrow(ATL_ACD))
ATL_ACD = bind_cols(ATL_ACD, 
          label1) %>% 
  rename(comparison = ...8)

ARC_ACD = read_tsv('ARC_ACD_Fst_100kb.txt')
label2 = rep('ARC-ACD', nrow(ARC_ACD))
ARC_ACD = bind_cols(ARC_ACD, 
                    label2) %>% 
  rename(comparison = ...8)

ARC_ATL = read_tsv('ARC_ATL_Fst_100kb.txt')
label3 = rep('ARC-ATL', nrow(ARC_ATL))
ARC_ATL = bind_cols(ARC_ATL, 
                    label3) %>% 
  rename(comparison = ...8)


glacial_lin_compare = ggplot() +
  geom_violin(data = ATL_ACD, 
              aes(y = FST_zero_mean, 
                  x = comparison), 
              fill = '#ff006e', 
              col = '#ff006e')+
  geom_point(data = ATL_ACD, 
             aes(y = 0.175105, 
                 x = comparison), 
             size = 3)+
  geom_violin(data = ARC_ACD, 
              aes(y = FST_zero_mean, 
                  x = comparison), 
              fill = '#fb8500', 
              col = '#fb8500')+
  geom_point(data = ARC_ACD, 
             aes(y = 0.397116, 
                 x = comparison), 
             size = 3)+
  geom_violin(data = ARC_ATL, 
              aes(y = FST_zero_mean, 
                  x = comparison), 
              fill = '#219ebc', 
              col = '#219ebc')+
  geom_point(data = ARC_ATL, 
             aes(y = 0.430073, 
                 x = comparison), 
             size = 3)+
  labs(y = 'Mean Fst per 100Kb')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 12))


ggsave('Glacial_lineage_Fst_100Kb.tiff', 
       plot = glacial_lin_compare, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)
