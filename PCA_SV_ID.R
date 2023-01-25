##############################
## Putative structural variant reanalysis
##
## Matt Brachmann (PhDMattyB)
##
## 2023-01-10
##
##############################

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
       path = '~/Bradbury_Postdoc/AC_SV_MS_Data/Figures/', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# detect sv per chr -------------------------------------------------------

setwd('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/Per_Chr_NewChrSet/')
identifiers = read_csv('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/Charr_Poly_All_Fixed_notbed_2.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID')

##
# LG1 SV detect -----------------------------------------------------------

LG1 = read.pcadapt('AC_New_CHRSET_1.bed', 
                     type = 'bed')

pca_LG1 = pcadapt(LG1, 
                    K = 10, 
                    method = 'mahalanobis', 
                    min.maf = 0.01)

plot(pca_LG1, 
     option = 'screeplot')
# 
# plot(pca_LG1, 
#      option = 'scores')

pca_LG1$singular.values
sum(pca_LG1$singular.values)

(sqrt(0.5733557)/1.968103)*100
## 38.47%
(sqrt(0.2538209)/1.968103)*100
## 25.60%
(sqrt(0.2335121)/1.968103)*100
## 24.55%

pca_LG1_scores = as_tibble(pca_LG1$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG1_scores.csv')

identifiers = read_csv('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('~/Bradbury_Postdoc/AC_SV_MS_Data/Pcadapt/Charr_Poly_All_Fixed_notbed_1.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG1_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'Chr 1')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr1.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG2 SV detect -----------------------------------------------------------

Chr2 = read.pcadapt('AC_New_CHRSET_2.bed', 
                     type = 'bed')

pca_Chr2 = pcadapt(Chr2, 
                    K = 10, 
                    method = 'mahalanobis', 
                    min.maf = 0.01)

plot(pca_Chr2, 
     option = 'screeplot')
# 
# plot(pca_Chr2, 
#      option = 'scores')
# summary(pca_Chr2)
# 
# pca_Chr2$singular.values
# sum(pca_Chr2$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_Chr2_scores = as_tibble(pca_Chr2$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr2_scores.csv')

meta_data = bind_cols(meta_data, 
            pca_Chr2_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'Chr2')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr2.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr3 SV detect -----------------------------------------------------------

Chr3 = read.pcadapt('AC_New_CHRSET_3.bed', 
                   type = 'bed')

pca_Chr3 = pcadapt(Chr3, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_Chr3, 
     option = 'screeplot')
# 
# plot(pca_Chr3, 
#      option = 'scores')
# summary(pca_Chr3)
# 
# pca_Chr3$singular.values
# sum(pca_Chr3$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_Chr3_scores = as_tibble(pca_Chr3$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr3_scores.csv')

chr3_meta =  bind_cols(meta_data, 
            pca_Chr3_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr3_meta %>%  
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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'Chr3')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr3.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr4 SV detect -----------------------------------------------------------

Chr4 = read.pcadapt('AC_New_CHRSET_4.bed', 
                   type = 'bed')

pca_Chr4 = pcadapt(Chr4, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_Chr4, 
     option = 'screeplot')
# 
# plot(pca_Chr4, 
#      option = 'scores')
# summary(pca_Chr4)
# 
# pca_Chr4$singular.values
# sum(pca_Chr4$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_Chr4_scores = as_tibble(pca_Chr4$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr4_scores.csv')

chr4_meta =  bind_cols(meta_data, 
            pca_Chr4_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr4_meta %>%  
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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'Chr4')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr4.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr5 SV detect ------------------------------------------------------

Chr5 = read.pcadapt('AC_New_CHRSET_5.bed', 
                   type = 'bed')

pca_chr5 = pcadapt(Chr5, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

# plot(pca_LG4q, 
#      option = 'screeplot')
# 
# plot(pca_LG4q, 
#      option = 'scores')
# summary(pca_LG4q)
# 
# pca_LG4q$singular.values
# sum(pca_LG4q$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_chr5_scores = as_tibble(pca_chr5$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_chr5_scores.csv')


chr5_meta = bind_cols(meta_data, 
            pca_chr5_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr5_meta %>%  
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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'Chr5')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr5.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr6 SV detect --------------------------------------------------------

Chr6 = read.pcadapt('AC_New_CHRSET_6.bed', 
                   type = 'bed')

pca_Chr6 = pcadapt(Chr6, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr6_scores = as_tibble(pca_Chr6$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr6_scores.csv')

chr6_meta = bind_cols(meta_data, 
            pca_Chr6_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr6_meta %>%  
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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'Chr6')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr6.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr7 SV detect -----------------------------------------------------------

Chr7 = read.pcadapt('AC_New_CHRSET_7.bed', 
                   type = 'bed')

pca_Chr7 = pcadapt(Chr7, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr7_scores = as_tibble(pca_Chr7$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr7_scores.csv')

chr7_meta = bind_cols(meta_data, 
            pca_Chr7_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr7_meta %>%  
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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'Chr7')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr7.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG6.1 sv detect ---------------------------------------------------------

LG6.1 = read.pcadapt('Charr_Poly_All_Fixed_8.bed', 
                   type = 'bed')

pca_LG6.1 = pcadapt(LG6.1, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG6.1, 
     option = 'screeplot')
# 
# plot(pca_LG6.1, 
#      option = 'scores')
# summary(pca_LG6.1)
# 
# pca_LG6.1$singular.values
# sum(pca_LG6.1$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG6.1_scores = as_tibble(pca_LG6.1$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG6.1_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_8.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG6.1_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG6.1')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG6.1.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG6.2 SV detect ------------------------------------------------

# map = read_tsv('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.map')

LG6.2 = read.pcadapt('Charr_Poly_All_Fixed_9.bed', 
                     type = 'bed')

pca_LG6.2 = pcadapt(LG6.2, 
                      K = 10, 
                      method = 'mahalanobis', 
                      min.maf = 0.01)

plot(pca_LG6.2, 
     option = 'screeplot')
# 
# plot(pca_LG6.2, 
#      option = 'scores')
# summary(pca_LG6.2)
# 
# pca_LG6.2$singular.values
# sum(pca_LG6.2$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG6.2_scores = as_tibble(pca_LG6.2$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG6.2_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_9.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG6.2_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG6.2')


ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG6.2.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)



# LG7 sv detect -----------------------------------------------------------

LG7 = read.pcadapt('Charr_Poly_All_Fixed_10.bed', 
                   type = 'bed')

pca_LG7 = pcadapt(LG7, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG7, 
     option = 'screeplot')
# 
# plot(pca_LG7, 
#      option = 'scores')
# summary(pca_LG7)
# 
# pca_LG7$singular.values
# sum(pca_LG7$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG7_scores = as_tibble(pca_LG7$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG7_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_10.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG7_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG7')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG7.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# AC08 or LG11 SV detect --------------------------------------------------

## Need to pull out chromosome LG11 from the full data set
## Then make the bed file(s) for pcadapt
LG11 = read.pcadapt('Charr_Poly_All_Fixed_11.bed', 
                     type = 'bed')

pca_LG11 = pcadapt(LG11, 
                    K = 10, 
                    method = 'mahalanobis', 
                    min.maf = 0.01)

plot(pca_LG11, 
     option = 'screeplot')
# 
# plot(pca_LG11, 
#      option = 'scores')
summary(pca_LG11)

pca_LG11$singular.values
sum(pca_LG11$singular.values)

(sqrt(0.5733557)/1.968103)*100
## 38.47%
(sqrt(0.2538209)/1.968103)*100
## 25.60%
(sqrt(0.2335121)/1.968103)*100
## 24.55%

pca_LG11_scores = as_tibble(pca_LG11$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG11_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_11.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG11_scores)

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
  labs(x = 'Principal component 1 ',
       y = 'Principal component 2 ', 
       col = 'Glacial lineages', 
       title = 'LG8')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG8.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)
# LG9 sv detect -----------------------------------------------------------

LG9 = read.pcadapt('Charr_Poly_All_Fixed_12.bed', 
                   type = 'bed')

pca_LG9 = pcadapt(LG9, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG9, 
     option = 'screeplot')
# 
# plot(pca_LG9, 
#      option = 'scores')
# summary(pca_LG9)
# 
# pca_LG9$singular.values
# sum(pca_LG9$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG9_scores = as_tibble(pca_LG9$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG9_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_12.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG9_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG9')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG9.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# LG10 sv detect ----------------------------------------------------------

LG10 = read.pcadapt('Charr_Poly_All_Fixed_13.bed', 
                   type = 'bed')

pca_LG10 = pcadapt(LG10, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG10, 
     option = 'screeplot')
# 
# plot(pca_LG10, 
#      option = 'scores')
# summary(pca_LG10)
# 
# pca_LG10$singular.values
# sum(pca_LG10$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG10_scores = as_tibble(pca_LG10$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG10_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_13.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG10_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG10')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG10.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)



# LG11 sv detect ----------------------------------------------------------

LG11 = read.pcadapt('Charr_Poly_All_Fixed_14.bed', 
                   type = 'bed')

pca_LG11 = pcadapt(LG11, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG11, 
     option = 'screeplot')
# 
# plot(pca_LG11, 
#      option = 'scores')
# summary(pca_LG11)
# 
# pca_LG11$singular.values
# sum(pca_LG11$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG11_scores = as_tibble(pca_LG11$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG11_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_14.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG11_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG11')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG11.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# LG12 --------------------------------------------------------------------

LG12 = read.pcadapt('Charr_Poly_All_Fixed_15.bed', 
                   type = 'bed')

pca_LG12 = pcadapt(LG12, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG12, 
     option = 'screeplot')
# 
# plot(pca_LG12, 
#      option = 'scores')
# summary(pca_LG12)
# 
# pca_LG12$singular.values
# sum(pca_LG12$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG12_scores = as_tibble(pca_LG12$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG12_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_15.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG12_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG12')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG12.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# LG13 sv detect ----------------------------------------------------------

LG13 = read.pcadapt('Charr_Poly_All_Fixed_16.bed', 
                   type = 'bed')

pca_LG13 = pcadapt(LG13, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG13, 
     option = 'screeplot')
# 
# plot(pca_LG13, 
#      option = 'scores')
# summary(pca_LG13)
# 
# pca_LG13$singular.values
# sum(pca_LG13$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG13_scores = as_tibble(pca_LG13$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG13_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_16.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG13_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG13')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG13.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG14 sv detect ----------------------------------------------------------

LG14 = read.pcadapt('Charr_Poly_All_Fixed_17.bed', 
                   type = 'bed')

pca_LG14 = pcadapt(LG14, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG14, 
     option = 'screeplot')
# 
# plot(pca_LG14, 
#      option = 'scores')
# summary(pca_LG14)
# 
# pca_LG14$singular.values
# sum(pca_LG14$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG14_scores = as_tibble(pca_LG14$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG14_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_17.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG14_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG14')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG14.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# LG15 --------------------------------------------------------------------

LG15 = read.pcadapt('Charr_Poly_All_Fixed_18.bed', 
                   type = 'bed')

pca_LG15 = pcadapt(LG15, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG15, 
     option = 'screeplot')
# 
# plot(pca_LG15, 
#      option = 'scores')
# summary(pca_LG15)
# 
# pca_LG15$singular.values
# sum(pca_LG15$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG15_scores = as_tibble(pca_LG15$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG15_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_18.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG15_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG15')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG15.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# LG16 --------------------------------------------------------------------

LG16 = read.pcadapt('Charr_Poly_All_Fixed_19.bed', 
                   type = 'bed')

pca_LG16 = pcadapt(LG16, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG16, 
     option = 'screeplot')
# 
# plot(pca_LG16, 
#      option = 'scores')
# summary(pca_LG16)
# 
# pca_LG16$singular.values
# sum(pca_LG16$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG16_scores = as_tibble(pca_LG16$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG16_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_19.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG16_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG16')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG16.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG17 sv detect ----------------------------------------------------------

LG17 = read.pcadapt('Charr_Poly_All_Fixed_20.bed', 
                   type = 'bed')

pca_LG17 = pcadapt(LG17, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG17, 
     option = 'screeplot')
# 
# plot(pca_LG17, 
#      option = 'scores')
# summary(pca_LG17)
# 
# pca_LG17$singular.values
# sum(pca_LG17$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG17_scores = as_tibble(pca_LG17$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG17_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_20.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG17_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG17')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG17.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG18 sv detect ----------------------------------------------------------

LG18 = read.pcadapt('Charr_Poly_All_Fixed_21.bed', 
                   type = 'bed')

pca_LG18 = pcadapt(LG18, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG18, 
     option = 'screeplot')
# 
# plot(pca_LG18, 
#      option = 'scores')
# summary(pca_LG18)
# 
# pca_LG18$singular.values
# sum(pca_LG18$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG18_scores = as_tibble(pca_LG18$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG18_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_21.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG18_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG18')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG18.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# LG19 sv detect ----------------------------------------------------------

LG19 = read.pcadapt('Charr_Poly_All_Fixed_22.bed', 
                   type = 'bed')

pca_LG19 = pcadapt(LG19, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG19, 
     option = 'screeplot')
# 
# plot(pca_LG19, 
#      option = 'scores')
# summary(pca_LG19)
# 
# pca_LG19$singular.values
# sum(pca_LG19$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG19_scores = as_tibble(pca_LG19$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG19_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_22.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG19_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG19')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG19.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG20 sv detect --------------------------------------------------------------------
LG20 = read.pcadapt('Charr_Poly_All_Fixed_23.bed', 
                   type = 'bed')

pca_LG20 = pcadapt(LG20, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG20, 
     option = 'screeplot')
# 
# plot(pca_LG20, 
#      option = 'scores')
# summary(pca_LG20)
# 
# pca_LG20$singular.values
# sum(pca_LG20$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG20_scores = as_tibble(pca_LG20$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG20_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_23.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG20_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG20')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG20.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG21 sv detect ----------------------------------------------------------

LG21 = read.pcadapt('Charr_Poly_All_Fixed_24.bed', 
                   type = 'bed')

pca_LG21 = pcadapt(LG21, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG21, 
     option = 'screeplot')
# 
# plot(pca_LG21, 
#      option = 'scores')
# summary(pca_LG21)
# 
# pca_LG21$singular.values
# sum(pca_LG21$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG21_scores = as_tibble(pca_LG21$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG21_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_24.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG21_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG21')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG21.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG22 --------------------------------------------------------------------

LG22 = read.pcadapt('Charr_Poly_All_Fixed_25.bed', 
                   type = 'bed')

pca_LG22 = pcadapt(LG22, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG22, 
     option = 'screeplot')
# 
# plot(pca_LG22, 
#      option = 'scores')
# summary(pca_LG22)
# 
# pca_LG22$singular.values
# sum(pca_LG22$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG22_scores = as_tibble(pca_LG22$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG22_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_25.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG22_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG22')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG22.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG23 sv detection -------------------------------------------------------

LG23 = read.pcadapt('Charr_Poly_All_Fixed_26.bed', 
                   type = 'bed')

pca_LG23 = pcadapt(LG23, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG23, 
     option = 'screeplot')
# 
# plot(pca_LG23, 
#      option = 'scores')
# summary(pca_LG23)
# 
# pca_LG23$singular.values
# sum(pca_LG23$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG23_scores = as_tibble(pca_LG23$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG23_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_26.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG23_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG23')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG23.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# LG24 sv detection -------------------------------------------------------

LG24 = read.pcadapt('Charr_Poly_All_Fixed_27.bed', 
                   type = 'bed')

pca_LG24 = pcadapt(LG24, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG24, 
     option = 'screeplot')
# 
# plot(pca_LG24, 
#      option = 'scores')
# summary(pca_LG24)
# 
# pca_LG24$singular.values
# sum(pca_LG24$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG24_scores = as_tibble(pca_LG24$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG24_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_27.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG24_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG24')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG24.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG25 sv detection -------------------------------------------------------

LG25 = read.pcadapt('Charr_Poly_All_Fixed_28.bed', 
                   type = 'bed')

pca_LG25 = pcadapt(LG25, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG25, 
     option = 'screeplot')
# 
# plot(pca_LG25, 
#      option = 'scores')
# summary(pca_LG25)
# 
# pca_LG25$singular.values
# sum(pca_LG25$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG25_scores = as_tibble(pca_LG25$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG25_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_28.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG25_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG25')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG25.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# LG26 sv detect ----------------------------------------------------------

LG26 = read.pcadapt('Charr_Poly_All_Fixed_29.bed', 
                   type = 'bed')

pca_LG26 = pcadapt(LG26, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG26, 
     option = 'screeplot')
# 
# plot(pca_LG26, 
#      option = 'scores')
# summary(pca_LG26)
# 
# pca_LG26$singular.values
# sum(pca_LG26$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG26_scores = as_tibble(pca_LG26$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG26_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_29.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG26_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG26')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG26.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# LG27 sv detect ----------------------------------------------------------

LG27 = read.pcadapt('Charr_Poly_All_Fixed_30.bed', 
                   type = 'bed')

pca_LG27 = pcadapt(LG27, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG27, 
     option = 'screeplot')
# 
# plot(pca_LG27, 
#      option = 'scores')
# summary(pca_LG27)
# 
# pca_LG27$singular.values
# sum(pca_LG27$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG27_scores = as_tibble(pca_LG27$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG27_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_30.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG27_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG27')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG27.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG28 sv detect ----------------------------------------------------------

LG28 = read.pcadapt('Charr_Poly_All_Fixed_31.bed', 
                   type = 'bed')

pca_LG28 = pcadapt(LG28, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG28, 
     option = 'screeplot')
# 
# plot(pca_LG28, 
#      option = 'scores')
# summary(pca_LG28)
# 
# pca_LG28$singular.values
# sum(pca_LG28$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG28_scores = as_tibble(pca_LG28$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG28_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_31.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG28_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG28')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG28.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG 29 sv detect -------------------------------------------------------------------

LG29 = read.pcadapt('Charr_Poly_All_Fixed_32.bed', 
                   type = 'bed')

pca_LG29 = pcadapt(LG29, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG29, 
     option = 'screeplot')
# 
# plot(pca_LG29, 
#      option = 'scores')
# summary(pca_LG29)
# 
# pca_LG29$singular.values
# sum(pca_LG29$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG29_scores = as_tibble(pca_LG29$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG29_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_32.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG29_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG29')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG29.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG30 sv detect ----------------------------------------------------------

LG30 = read.pcadapt('Charr_Poly_All_Fixed_33.bed', 
                   type = 'bed')

pca_LG30 = pcadapt(LG30, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG30, 
     option = 'screeplot')
# 
# plot(pca_LG30, 
#      option = 'scores')
# summary(pca_LG30)
# 
# pca_LG30$singular.values
# sum(pca_LG30$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG30_scores = as_tibble(pca_LG30$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG30_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_33.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG30_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG30')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG30.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG31 sv detect ----------------------------------------------------------

LG31 = read.pcadapt('Charr_Poly_All_Fixed_34.bed', 
                   type = 'bed')

pca_LG31 = pcadapt(LG31, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG31, 
     option = 'screeplot')
# 
# plot(pca_LG31, 
#      option = 'scores')
# summary(pca_LG31)
# 
# pca_LG31$singular.values
# sum(pca_LG31$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG31_scores = as_tibble(pca_LG31$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG31_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_34.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG31_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG31')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG31.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG32 sv detection -------------------------------------------------------

LG32 = read.pcadapt('Charr_Poly_All_Fixed_35.bed', 
                   type = 'bed')

pca_LG32 = pcadapt(LG32, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG32, 
     option = 'screeplot')
# 
# plot(pca_LG32, 
#      option = 'scores')
# summary(pca_LG32)
# 
# pca_LG32$singular.values
# sum(pca_LG32$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG32_scores = as_tibble(pca_LG32$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG32_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_35.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG32_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG32')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG32.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG33 sv detect ----------------------------------------------------------

LG33 = read.pcadapt('Charr_Poly_All_Fixed_36.bed', 
                   type = 'bed')

pca_LG33 = pcadapt(LG33, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG33, 
     option = 'screeplot')
# 
# plot(pca_LG33, 
#      option = 'scores')
# summary(pca_LG33)
# 
# pca_LG33$singular.values
# sum(pca_LG33$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG33_scores = as_tibble(pca_LG33$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG33_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_36.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG33_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG33')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG33.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG34 sv detect ----------------------------------------------------------

LG34 = read.pcadapt('Charr_Poly_All_Fixed_37.bed', 
                   type = 'bed')

pca_LG34 = pcadapt(LG34, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG34, 
     option = 'screeplot')
# 
# plot(pca_LG34, 
#      option = 'scores')
# summary(pca_LG34)
# 
# pca_LG34$singular.values
# sum(pca_LG34$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG34_scores = as_tibble(pca_LG34$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG34_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_37.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG34_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG34')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG34.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG35 sv detect ----------------------------------------------------------

LG35 = read.pcadapt('Charr_Poly_All_Fixed_38.bed', 
                   type = 'bed')

pca_LG35 = pcadapt(LG35, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG35, 
     option = 'screeplot')
# 
# plot(pca_LG35, 
#      option = 'scores')
# summary(pca_LG35)
# 
# pca_LG35$singular.values
# sum(pca_LG35$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG35_scores = as_tibble(pca_LG35$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG35_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_38.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG35_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG35')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG35.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG36 sv detect ----------------------------------------------------------

LG36 = read.pcadapt('Charr_Poly_All_Fixed_39.bed', 
                   type = 'bed')

pca_LG36 = pcadapt(LG36, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG36, 
     option = 'screeplot')
# 
# plot(pca_LG36, 
#      option = 'scores')
# summary(pca_LG36)
# 
# pca_LG36$singular.values
# sum(pca_LG36$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG36_scores = as_tibble(pca_LG36$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG36_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_39.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG36_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG36')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG36.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG37  sv detection--------------------------------------------------------------------

LG37 = read.pcadapt('Charr_Poly_All_Fixed_40.bed', 
                   type = 'bed')

pca_LG37 = pcadapt(LG37, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG37, 
     option = 'screeplot')
# 
# plot(pca_LG37, 
#      option = 'scores')
# summary(pca_LG37)
# 
# pca_LG37$singular.values
# sum(pca_LG37$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG37_scores = as_tibble(pca_LG37$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG37_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_40.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG37_scores)

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
  labs(x = 'Principal component 1',
       y = 'Principal component 2', 
       col = 'Glacial lineages', 
       title = 'LG37')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG37.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)






