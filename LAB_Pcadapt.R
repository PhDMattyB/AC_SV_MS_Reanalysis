

setwd('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/')

library(data.table)
library(tidyverse)
library(pcadapt)
library(devtools)
# install_github('isglobal-brge/invClust')
library(invClust)
library(patchwork)

theme_set(theme_bw())



# All populations pca -----------------------------------------------------

allpops = read.pcadapt('Charr_Poly_All_Fixed_coords_maf05_geno95_envmatch_Lab.bed', 
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

(sqrt(0.26759009)/1.141238)*100
## 45.32%
(sqrt(0.13393490)/1.141238)*100
## 32.07%
(sqrt(0.11862102)/1.141238)*100
## 30.18%

pca_allpops_scores = as_tibble(pca_allpops$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('LAB_Charr_PCA_allpops_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Lab_Charr.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_allpops_scores)


# glacial_cols = c('#F23545',
#                           '#4E458C',
#                           '#4E94BF', 
#                           '#F29F05')
                          
meta_data %>% 
  group_by(FID)

meta_data  %>% 
  ggplot(aes(x = PC1, y = PC2))+
  # geom_point(aes(col = Location), 
  #            size = 2)+
  geom_point(aes(col = Lat), 
             size = 2)+
  # geom_point(aes(col = Population), 
  #            size = 2)+
  # scale_color_manual(values = cols)+
  # scale_color_manual(values = glacial_cols)+
  scale_color_viridis(option = 'magma')+
  scale_fill_viridis(option = 'magma')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(size = 15, 
                                  hjust = 0), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.position = 'none') +
  labs(x = 'Principal component 1 (45.3%)',
       y = 'Principal component 2 (32.1%)', 
       col = 'Glacial lineages')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages.tiff',
       path = '~/Bradbury_Postdoc/AC_SV_MS_Data/Figures/', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# detect sv per chr -------------------------------------------------------

setwd('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/LAB_PCA_per_Chr/')
identifiers = read_csv('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/Lab_Charr.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID')

# Chr13 sv detect ----------------------------------------------------------

Chr13 = read.pcadapt('Lab_Charr_13.bed', 
                     type = 'bed')

pca_Chr13 = pcadapt(Chr13, 
                    K = 10, 
                    method = 'mahalanobis', 
                    min.maf = 0.01)

pca_Chr13_scores = as_tibble(pca_Chr13$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Lab_Charr_PCA_Chr13_scores.csv')

chr13_meta = bind_cols(meta_data, 
                       pca_Chr13_scores)


chr13_meta %>%  
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(aes(col = Lat), 
             size = 2)+
  # geom_point(aes(col = Location), 
  #            size = 2)+
  # geom_point(aes(col = Glacial_lin), 
  #            size = 2)+
  # geom_point(aes(col = Population), 
  #            size = 2)+
  # scale_color_manual(values = cols)+
  scale_colour_viridis(option = 'magma')+
  # scale_color_manual(values = glacial_cols)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(size = 15, 
                                  hjust = 0), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.position = 'none') +
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'Chr13')

ggsave(file = 'LAB_PCAdapt__Chr13.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)
