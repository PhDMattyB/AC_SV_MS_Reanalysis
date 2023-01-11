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
library(patchwork)

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

meta_data = read_table2('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.ped', 
              col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_allpops_scores)


glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
meta_data %>%  
  ggplot(aes(x = PC1, y = PC2))+
  # geom_point(aes(col = Location), 
  #            size = 2)+
  geom_point(aes(col = Glacial_lin), 
             size = 2)+
  # geom_point(aes(col = Population), 
  #            size = 2)+
  # scale_color_manual(values = cols)+
  scale_color_manual(values = glacial_cols)+
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
  labs(x = 'Principal component 1 (43.7%)',
       y = 'Principal component 2 (29.6%)', 
       col = 'Glacial lineages')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# AC08 or LG11 SV detect --------------------------------------------------

## Need to pull out chromosome LG11 from the full data set


