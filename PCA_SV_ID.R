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

# LG1 SV detect -----------------------------------------------------------

LG1 = read.pcadapt('Charr_Poly_All_Fixed_1.bed', 
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

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_1.ped', 
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
       title = 'LG1')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG1.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG2 SV detect -----------------------------------------------------------

LG2 = read.pcadapt('Charr_Poly_All_Fixed_2.bed', 
                     type = 'bed')

pca_LG2 = pcadapt(LG2, 
                    K = 10, 
                    method = 'mahalanobis', 
                    min.maf = 0.01)

plot(pca_LG2, 
     option = 'screeplot')
# 
# plot(pca_LG2, 
#      option = 'scores')
# summary(pca_LG2)
# 
# pca_LG2$singular.values
# sum(pca_LG2$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG2_scores = as_tibble(pca_LG2$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG2_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_2.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG2_scores)

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
       title = 'LG2')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG2.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG3 SV detect -----------------------------------------------------------

LG3 = read.pcadapt('Charr_Poly_All_Fixed_3.bed', 
                   type = 'bed')

pca_LG3 = pcadapt(LG3, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG3, 
     option = 'screeplot')
# 
# plot(pca_LG3, 
#      option = 'scores')
# summary(pca_LG3)
# 
# pca_LG3$singular.values
# sum(pca_LG3$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG3_scores = as_tibble(pca_LG3$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG3_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_3.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG3_scores)

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
       title = 'LG3')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG3.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG4 SV detect -----------------------------------------------------------

LG4 = read.pcadapt('Charr_Poly_All_Fixed_4.bed', 
                   type = 'bed')

pca_LG4 = pcadapt(LG4, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG4, 
     option = 'screeplot')
# 
# plot(pca_LG4, 
#      option = 'scores')
# summary(pca_LG4)
# 
# pca_LG4$singular.values
# sum(pca_LG4$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG4_scores = as_tibble(pca_LG4$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG4_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_4.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG4_scores)

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
       title = 'LG4p')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG4p.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG4q1:29 SV detect ------------------------------------------------------

LG4q = read.pcadapt('Charr_Poly_All_Fixed_5.bed', 
                   type = 'bed')

pca_LG4q = pcadapt(LG4q, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG4q, 
     option = 'screeplot')
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

pca_LG4q_scores = as_tibble(pca_LG4q$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG4q_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_5.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG4q_scores)

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
       title = 'LG4q1:29')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG4q1:29.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG4q.2 SV detect --------------------------------------------------------

LG4q.2 = read.pcadapt('Charr_Poly_All_Fixed_6.bed', 
                   type = 'bed')

pca_LG4q.2 = pcadapt(LG4q.2, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG4q.2, 
     option = 'screeplot')
# 
# plot(pca_LG4q.2, 
#      option = 'scores')
# summary(pca_LG4q.2)
# 
# pca_LG4q.2$singular.values
# sum(pca_LG4q.2$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG4q.2_scores = as_tibble(pca_LG4q.2$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG4q.2_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_6.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG4q.2_scores)

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
       title = 'LG4q.2')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG4q.2.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


# LG5 SV detect -----------------------------------------------------------

LG5 = read.pcadapt('Charr_Poly_All_Fixed_7.bed', 
                   type = 'bed')

pca_LG5 = pcadapt(LG5, 
                  K = 10, 
                  method = 'mahalanobis', 
                  min.maf = 0.01)

plot(pca_LG5, 
     option = 'screeplot')
# 
# plot(pca_LG5, 
#      option = 'scores')
# summary(pca_LG5)
# 
# pca_LG5$singular.values
# sum(pca_LG5$singular.values)
# 
# (sqrt(0.5733557)/1.968103)*100
# ## 38.47%
# (sqrt(0.2538209)/1.968103)*100
# ## 25.60%
# (sqrt(0.2335121)/1.968103)*100
# ## 24.55%

pca_LG5_scores = as_tibble(pca_LG5$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3) %>%
  dplyr::select(PC1,
                PC2,
                PC3) %>% 
  write_csv('Charr_PCA_LG5_scores.csv')

identifiers = read_csv('ggtree_labels.csv') %>% 
  rename(FID = Population)

meta_data = read_table2('Charr_Poly_All_Fixed_notbed_7.ped', 
                        col_names = F) %>% 
  dplyr::select(X1:X2) %>% 
  rename(FID = X1, 
         IID = X2) %>% 
  left_join(., 
            identifiers, 
            by = 'FID') %>% 
  bind_cols(., 
            pca_LG5_scores)

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
       title = 'LG5')

ggsave(file = 'PCAdapt_all_pops_k4_Glacial_lineages_LG5.tiff', 
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





