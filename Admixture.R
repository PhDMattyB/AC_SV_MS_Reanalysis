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
library(reshape2)

qvalues = read_table('Charr_Poly_All_Fixed_coords_maf05_geno95.4.Q', 
                   col_names = c('Q1', 
                                 'Q2', 
                                 'Q3', 
                                 'Q4'))

ids = read_table2('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.ped', 
            col_names = F) %>% 
  dplyr::select(1:2)

Clean_data = bind_cols(ids, 
                       qvalues) %>% 
  rename(Population = X1, 
         Individual = X2)

latlong = read_csv('SampleSites_Coords_1June2020.csv')


Clean_data = inner_join(Clean_data, 
           latlong) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin, 
                Lat, 
                Long, 
                Q1, 
                Q2, 
                Q3, 
                Q4)

# Graph admixture ---------------------------------------------------------


melted_dwata = melt(Clean_data, 
                    id.vars = c('Population', 
                                'Individual', 
                                'Glacial_lin', 
                                'Lat',
                                'Long')) %>% 
  as_tibble() 

## bw is the king of plots
theme_set(theme_bw())

## need a colour palette
test_col = c( '#4E9EBF',
              '#4E458C',
              '#F29F05',
              '#F23545')


admixture = ggplot(data = melted_dwata, 
                      aes(x = reorder(Individual, 
                                      Lat),
                          y = value, 
                          fill = variable))+
  geom_bar(stat = "identity", 
           width = 1)+
  scale_fill_manual(values = test_col)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        ## can add xaxis labels if needed
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 8,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

admixture

ggsave('admixture_k4_glacial_lineages.tiff',
       plot = admixture, 
       dpi = 'retina', 
       units = 'cm', 
       width = 35, 
       height = 10)


# admixture bioclim outliers ----------------------------------------------

bioclim_admix = read_table2('bioclim_outlier_snps.3.Q', 
                         col_names = c('Q1', 
                                       'Q2', 
                                       'Q3'))

bioclim_data = bind_cols(ids, 
                         bioclim_admix) %>% 
  rename(Population = X1, 
         Individual = X2)

bioclim_data = inner_join(bioclim_data, 
                         latlong) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin, 
                Lat, 
                Long, 
                Q1, 
                Q2, 
                Q3)


bioclim_melted = melt(bioclim_data, 
                    id.vars = c('Population', 
                                'Individual', 
                                'Glacial_lin', 
                                'Lat',
                                'Long')) %>% 
  as_tibble() 

## bw is the king of plots
theme_set(theme_bw())

## need a colour palette
# test_col = c( '#4E458C',
#               '#4E9EBF',
#               '#F23545',
#               '#F29F05')

test_col = c('#4E9EBF',
             '#F23545',
             '#4E458C')

bioclim_melted %>% View()

bioclim_admixture = ggplot(data = bioclim_melted, 
                   aes(x = reorder(Individual, 
                                   Lat),
                       y = value, 
                       fill = variable))+
  geom_bar(stat = "identity", 
           width = 1)+
  scale_fill_manual(values = test_col)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        ## can add xaxis labels if needed
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 8,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

bioclim_admixture

ggsave('admixture_k4_bioclim_outliers.tiff',
       plot = bioclim_admixture, 
       dpi = 'retina', 
       units = 'cm', 
       width = 35, 
       height = 10)

# admixture smd outleirs --------------------------------------------------

sdm_admix = read_table2('sdm_outlier_snps.3.Q', 
                            col_names = c('Q1', 
                                          'Q2', 
                                          'Q3'))

sdm_data = bind_cols(ids, 
                         sdm_admix) %>% 
  rename(Population = X1, 
         Individual = X2)

sdm_data = inner_join(sdm_data, 
                          latlong) %>% 
  dplyr::select(Population, 
                Individual, 
                Glacial_lin, 
                Lat, 
                Long, 
                Q1, 
                Q2, 
                Q3)


sdm_melted = melt(sdm_data, 
                      id.vars = c('Population', 
                                  'Individual', 
                                  'Glacial_lin', 
                                  'Lat',
                                  'Long')) %>% 
  as_tibble() 

## bw is the king of plots
theme_set(theme_bw())

## need a colour palette
# test_col = c( '#F23545',
#               '#F29F05',
#               '#4E9EBF',
#               '#4E458C')

test_col = c('#4E9EBF',
             '#F23545',
             '#4E458C')

sdm_melted %>% View()

sdm_admixture = ggplot(data = sdm_melted, 
                           aes(x = reorder(Individual, 
                                           Lat),
                               y = value, 
                               fill = variable))+
  geom_bar(stat = "identity", 
           width = 1)+
  scale_fill_manual(values = test_col)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        ## can add xaxis labels if needed
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 8,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

sdm_admixture

ggsave('admixture_k3_sdm_outliers.tiff',
       plot = sdm_admixture, 
       dpi = 'retina', 
       units = 'cm', 
       width = 35, 
       height = 10)


# snmf data cleaning ---------------------------------------------------------------

snmf_qvalues = read_csv('Charr_Poly_snmf_K4_qvalues.csv')

group = read_tsv('Charr_Poly_snmf_K4.txt')

snmf_data = bind_cols(group, 
                      snmf_qvalues)


ids = read_table2('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.ped', 
                  col_names = F) %>% 
  dplyr::select(1:2)

snmf_data = bind_cols(ids, 
                       snmf_data) %>% 
  rename(Population = X1, 
         Individual = X2)

latlong = read_csv('SampleSites_Coords_1June2020.csv')


snmf_data = inner_join(snmf_data, 
                        latlong) %>% 
  dplyr::select(Population, 
                Individual, 
                Genetic_group,
                Glacial_lin, 
                Lat, 
                Long, 
                Q1, 
                Q2, 
                Q3, 
                Q4)


# snmf plot ---------------------------------------------------------------

snmf_melted = melt(snmf_data, 
                    id.vars = c('Population', 
                                'Individual', 
                                'Genetic_group',
                                'Glacial_lin', 
                                'Lat',
                                'Long')) %>% 
  as_tibble() 
## bw is the king of plots
theme_set(theme_bw())

## need a colour palette
test_col = c( '#F23545',
              '#4E9EBF',
              '#F29F05',
              '#4E458C')


snmf_k4 = ggplot(data = snmf_melted, 
                   aes(x = reorder(Individual, 
                                   Lat),
                       y = value, 
                       fill = variable))+
  geom_bar(stat = "identity", 
           width = 1)+
  scale_fill_manual(values = test_col)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        ## can add xaxis labels if needed
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 8,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

snmf_k4

ggsave('snmf_k4_glacial_lineages.tiff',
       plot = snmf_k4, 
       dpi = 'retina', 
       units = 'cm', 
       width = 35, 
       height = 10)

