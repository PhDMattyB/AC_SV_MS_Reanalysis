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


# Chr8 sv detect ---------------------------------------------------------

Chr8 = read.pcadapt('AC_New_CHRSET_8.bed', 
                   type = 'bed')

pca_Chr8 = pcadapt(Chr8, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr8_scores = as_tibble(pca_Chr8$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr8_scores.csv')

chr8_meta = bind_cols(meta_data, 
            pca_Chr8_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr8_meta %>%  
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
       title = 'Chr8')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr8.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr9 SV detect ------------------------------------------------

# map = read_tsv('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.map')

Chr9 = read.pcadapt('AC_New_CHRSET_9.bed', 
                     type = 'bed')

pca_Chr9 = pcadapt(Chr9, 
                      K = 10, 
                      method = 'mahalanobis', 
                      min.maf = 0.01)

pca_Chr9_scores = as_tibble(pca_Chr9$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr9_scores.csv')

chr9_meta = bind_cols(meta_data, 
            pca_Chr9_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr9_meta %>%  
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
       title = 'Chr9')


ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr9.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)



# Chr10 sv detect -----------------------------------------------------------

Chr10 = read.pcadapt('AC_New_CHRSET_10.bed', 
                   type = 'bed')

pca_Chr10 = pcadapt(Chr10, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr10_scores = as_tibble(pca_Chr10$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr10_scores.csv')

chr10_meta = bind_cols(meta_data, 
            pca_Chr10_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr10_meta %>%  
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
       title = 'Chr10')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr10.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr11 SV detect --------------------------------------------------

## Need to pull out chromosome Chr11 from the full data set
## Then make the bed file(s) for pcadapt
Chr11 = read.pcadapt('AC_New_CHRSET_11.bed', 
                     type = 'bed')

pca_Chr11 = pcadapt(Chr11, 
                    K = 10, 
                    method = 'mahalanobis', 
                    min.maf = 0.01)

pca_Chr11_scores = as_tibble(pca_Chr11$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr11_scores.csv')

chr11_meta = bind_cols(meta_data, 
            pca_Chr11_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr11_meta %>%  
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
       title = 'Chr11')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr11.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)
# Chr12 sv detect -----------------------------------------------------------

Chr12 = read.pcadapt('AC_New_CHRSET_12.bed', 
                   type = 'bed')

pca_Chr12 = pcadapt(Chr12, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr12_scores = as_tibble(pca_Chr12$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr12_scores.csv')

chr12_meta = bind_cols(meta_data, 
            pca_Chr12_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr12_meta %>%  
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
       title = 'Chr12')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr12.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# Chr13 sv detect ----------------------------------------------------------

Chr13 = read.pcadapt('AC_New_CHRSET_13.bed', 
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
  write_csv('Charr_PCA_Chr13_scores.csv')

chr13_meta = bind_cols(meta_data, 
            pca_Chr13_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr13_meta %>%  
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
       title = 'Chr13')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr13.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)



# Chr14 sv detect ----------------------------------------------------------

Chr14 = read.pcadapt('AC_New_CHRSET_14.bed', 
                   type = 'bed')

pca_Chr14 = pcadapt(Chr14, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr14_scores = as_tibble(pca_Chr14$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr14_scores.csv')

chr14_meta = bind_cols(meta_data, 
            pca_Chr14_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr14_meta %>%  
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
       title = 'Chr14')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr14.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# Chr15 --------------------------------------------------------------------

Chr15 = read.pcadapt('AC_New_CHRSET_15.bed', 
                   type = 'bed')

pca_Chr15 = pcadapt(Chr15, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr15_scores = as_tibble(pca_Chr15$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr15_scores.csv')

chr15_meta = bind_cols(meta_data, 
            pca_Chr15_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr15_meta %>%  
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
       title = 'Chr15')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr15.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# Chr16 sv detect ----------------------------------------------------------

Chr16 = read.pcadapt('AC_New_CHRSET_16.bed', 
                   type = 'bed')

pca_Chr16 = pcadapt(Chr16, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr16_scores = as_tibble(pca_Chr16$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr16_scores.csv')

chr16_meta = bind_cols(meta_data, 
            pca_Chr16_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr16_meta %>%  
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
       title = 'Chr16')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr16.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr17 sv detect ----------------------------------------------------------

Chr17 = read.pcadapt('AC_New_CHRSET_17.bed', 
                   type = 'bed')

pca_Chr17 = pcadapt(Chr17, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr17_scores = as_tibble(pca_Chr17$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr17_scores.csv')

chr17_meta =  bind_cols(meta_data, 
            pca_Chr17_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr17_meta %>%  
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
       title = 'Chr17')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr17.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# Chr18 --------------------------------------------------------------------

Chr18 = read.pcadapt('AC_New_CHRSET_18.bed', 
                   type = 'bed')

pca_Chr18 = pcadapt(Chr18, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr18_scores = as_tibble(pca_Chr18$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr18_scores.csv')

chr18_meta =  bind_cols(meta_data, 
            pca_Chr18_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr18_meta %>%  
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
       title = 'Chr18')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr18.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# CHr19 --------------------------------------------------------------------

CHr19 = read.pcadapt('AC_New_CHRSET_19.bed', 
                   type = 'bed')

pca_CHr19 = pcadapt(CHr19, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_CHr19_scores = as_tibble(pca_CHr19$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_CHr19_scores.csv')

chr19_meta =  bind_cols(meta_data, 
            pca_CHr19_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr19_meta %>%  
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
       title = 'CHr19')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_CHr19.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr20 sv detect ----------------------------------------------------------

Chr20 = read.pcadapt('AC_New_CHRSET_20.bed', 
                   type = 'bed')

pca_Chr20 = pcadapt(Chr20, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr20_scores = as_tibble(pca_Chr20$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr20_scores.csv')

chr20_meta  = bind_cols(meta_data, 
            pca_Chr20_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr20_meta %>%  
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
       title = 'Chr20')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr20.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr21 sv detect ----------------------------------------------------------

Chr21 = read.pcadapt('AC_New_CHRSET_21.bed', 
                   type = 'bed')

pca_Chr21 = pcadapt(Chr21, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr21_scores = as_tibble(pca_Chr21$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr21_scores.csv')

chr21_meta  = bind_cols(meta_data, 
            pca_Chr21_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr21_meta %>%  
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
       title = 'Chr21')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr21.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# Chr22 sv detect ----------------------------------------------------------

Chr22 = read.pcadapt('AC_New_CHRSET_22.bed', 
                   type = 'bed')

pca_Chr22 = pcadapt(Chr22, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr22_scores = as_tibble(pca_Chr22$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr22_scores.csv')

chr22_meta =  bind_cols(meta_data, 
            pca_Chr22_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr22_meta %>%  
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
       title = 'Chr22')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr22.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr23 sv detect --------------------------------------------------------------------
Chr23 = read.pcadapt('AC_New_CHRSET_23.bed', 
                   type = 'bed')

pca_Chr23 = pcadapt(Chr23, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr23_scores = as_tibble(pca_Chr23$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr23_scores.csv')

chr23_meta =  bind_cols(meta_data, 
            pca_Chr23_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr23_meta %>%  
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
       title = 'Chr23')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr23.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr24 sv detect ----------------------------------------------------------

Chr24 = read.pcadapt('AC_New_CHRSET_24.bed', 
                   type = 'bed')

pca_Chr24 = pcadapt(Chr24, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr24_scores = as_tibble(pca_Chr24$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr24_scores.csv')

chr24_meta = bind_cols(meta_data, 
            pca_Chr24_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr24_meta %>%  
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
       title = 'Chr24')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr24.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr25 --------------------------------------------------------------------

Chr25 = read.pcadapt('AC_New_CHRSET_25.bed', 
                   type = 'bed')

pca_Chr25 = pcadapt(Chr25, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr25_scores = as_tibble(pca_Chr25$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr25_scores.csv')

chr25_meta =  bind_cols(meta_data, 
            pca_Chr25_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr25_meta %>%  
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
       title = 'Chr25')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr25.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr26 sv detection -------------------------------------------------------

Chr26 = read.pcadapt('AC_New_CHRSET_26.bed', 
                   type = 'bed')

pca_Chr26 = pcadapt(Chr26, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr26_scores = as_tibble(pca_Chr26$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr26_scores.csv')

chr26_meta =  bind_cols(meta_data, 
            pca_Chr26_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr26_meta %>%  
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
       title = 'Chr26')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr26.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# Chr27 sv detection -------------------------------------------------------

Chr27 = read.pcadapt('AC_New_CHRSET_27.bed', 
                   type = 'bed')

pca_Chr27 = pcadapt(Chr27, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr27_scores = as_tibble(pca_Chr27$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr27_scores.csv')

chr27_meta =  bind_cols(meta_data, 
            pca_Chr27_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr27_meta %>%  
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
       title = 'Chr27')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr27.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr28 sv detection -------------------------------------------------------

Chr28 = read.pcadapt('AC_New_CHRSET_28.bed', 
                   type = 'bed')

pca_Chr28 = pcadapt(Chr28, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr28_scores = as_tibble(pca_Chr28$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr28_scores.csv')

chr28_meta =  bind_cols(meta_data, 
            pca_Chr28_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr28_meta %>%  
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
       title = 'Chr28')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr28.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# Chr29 sv detect ----------------------------------------------------------

Chr29 = read.pcadapt('AC_New_CHRSET_29.bed', 
                   type = 'bed')

pca_Chr29 = pcadapt(Chr29, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr29_scores = as_tibble(pca_Chr29$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr29_scores.csv')

chr29_meta =  bind_cols(meta_data, 
            pca_Chr29_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr29_meta %>%  
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
       title = 'Chr29')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr29.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# Chr30 sv detect ----------------------------------------------------------

Chr30 = read.pcadapt('AC_New_CHRSET_30.bed', 
                   type = 'bed')

pca_Chr30 = pcadapt(Chr30, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr30_scores = as_tibble(pca_Chr30$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr30_scores.csv')

chr30_meta =  bind_cols(meta_data, 
            pca_Chr30_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr30_meta %>%  
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
       title = 'Chr30')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr30.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr31 sv detect ----------------------------------------------------------

Chr31 = read.pcadapt('AC_New_CHRSET_31.bed', 
                   type = 'bed')

pca_Chr31 = pcadapt(Chr31, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)


pca_Chr31_scores = as_tibble(pca_Chr31$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr31_scores.csv')

chr31_meta =  bind_cols(meta_data, 
            pca_Chr31_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr31_meta %>%  
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
       title = 'Chr31')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr31.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr32 sv detect -------------------------------------------------------------------

Chr32 = read.pcadapt('AC_New_CHRSET_32.bed', 
                   type = 'bed')

pca_Chr32 = pcadapt(Chr32, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr32_scores = as_tibble(pca_Chr32$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr32_scores.csv')

chr32_meta =  bind_cols(meta_data, 
            pca_Chr32_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr32_meta %>%  
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
       title = 'Chr32')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr32.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr33 sv detect ----------------------------------------------------------

Chr33 = read.pcadapt('AC_New_CHRSET_33.bed', 
                   type = 'bed')

pca_Chr33 = pcadapt(Chr33, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr33_scores = as_tibble(pca_Chr33$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr33_scores.csv')

chr33_meta =  bind_cols(meta_data, 
            pca_Chr33_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr33_meta %>%  
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
       title = 'Chr33')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr33.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr34 sv detect ----------------------------------------------------------

Chr34 = read.pcadapt('AC_New_CHRSET_34.bed', 
                   type = 'bed')

pca_Chr34 = pcadapt(Chr34, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr34_scores = as_tibble(pca_Chr34$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr34_scores.csv')

chr34_meta =  bind_cols(meta_data, 
            pca_Chr34_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr34_meta %>%  
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
       title = 'Chr34')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr34.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr35 sv detection -------------------------------------------------------

Chr35 = read.pcadapt('AC_New_CHRSET_35.bed', 
                   type = 'bed')

pca_Chr35 = pcadapt(Chr35, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

pca_Chr35_scores = as_tibble(pca_Chr35$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr35_scores.csv')

chr35_meta =  bind_cols(meta_data, 
            pca_Chr35_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr35_meta %>%  
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
       title = 'Chr35')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr35.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr36 sv detect ----------------------------------------------------------

Chr36 = read.pcadapt('AC_New_CHRSET_36.bed', 
                   type = 'bed')

pca_Chr36 = pcadapt(Chr36, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)


pca_Chr36_scores = as_tibble(pca_Chr36$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr36_scores.csv')

chr36_meta =  bind_cols(meta_data, 
            pca_Chr36_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr36_meta %>%  
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
       title = 'Chr36')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr36.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# Chr37 sv detect ----------------------------------------------------------

Chr37 = read.pcadapt('AC_New_CHRSET_37.bed', 
                   type = 'bed')

pca_Chr37 = pcadapt(Chr37, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)


pca_Chr37_scores = as_tibble(pca_Chr37$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_Chr37_scores.csv')

chr37_meta =  bind_cols(meta_data, 
            pca_Chr37_scores)

glacial_cols = c('#F23545',
                          '#4E458C',
                          '#4E94BF', 
                          '#F29F05')
                          
chr37_meta %>%  
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
       title = 'Chr37')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_Chr37.tiff', 
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






