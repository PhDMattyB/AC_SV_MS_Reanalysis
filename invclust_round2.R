##############################
## Invclust across each chr
##
## Matt Brachmann (PhDMattyB)
##
## 2023-01-14
##
##############################

setwd('~/Bradbury_Postdoc/AC_SV_MS_Data/invclust/')

library(patchwork)
library(invClust)
library(tidyverse)
library(snpStats)
# BiocManager::install("inveRsion")
# library(inveRsion)
# library(scoreInvHap)

theme_set(theme_bw())

## Not in function
'%!in%' <- function(x,y)!('%in%'(x,y))


# Chr5 --------------------------------------------------------------------
chr5 = read.plink(bed = 'AC_New_CHRSET_5.bed', 
                        bim = 'AC_New_CHRSET_5.bim', 
                        fam = 'AC_New_CHRSET_5.fam')

chr5_geno = chr5$genotypes
chr5_map = chr5$map 
identical(chr5_map[,2], 
          colnames(chr5_geno))

## SNP cluster 1
chr5_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr5 = data.frame(chr = 5, 
                              LBP = 5112629, 
                              RBP = 27224716, 
                              reg = 'inver1')

chr5_inver1 = invClust(roi = ROI_1_chr5, 
                             wh = 1, 
                             geno = chr5_geno, 
                             annot = chr5_map, 
                             dim = 2)

## SNP cluster 2
chr5_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr5 = data.frame(chr = 5, 
                              LBP = 27324207, 
                              RBP = 35666697, 
                              reg = 'inver2')

chr5_inver2 = invClust(roi = ROI_2_chr5, 
                             wh = 1, 
                             geno = chr5_geno, 
                             annot = chr5_map, 
                             dim = 2)

## SNP cluster 3
chr5_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr5 = data.frame(chr = 5, 
                              LBP = 35704209, 
                              RBP = 46880380, 
                              reg = 'inver3')

chr5_inver3 = invClust(roi = ROI_3_chr5, 
                             wh = 1, 
                             geno = chr5_geno, 
                             annot = chr5_map, 
                             dim = 2)

## SNP cluster 4
chr5_map %>% 
  slice(301:400) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_4_chr5 = data.frame(chr = 5,
                              LBP = 47107486,
                              RBP = 59083756,
                              reg = 'inver4')

chr5_inver4 = invClust(roi = ROI_4_chr5,
                             wh = 1,
                             geno = chr5_geno,
                             annot = chr5_map,
                             dim = 2)


## SNP cluster 5
chr5_map %>% 
  slice(401:nrow(chr5_map)) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_5_chr5 = data.frame(chr = 5,
                              LBP = 59107209,
                              RBP = 76481077,
                              reg = 'inver5')

chr5_inver5 = invClust(roi = ROI_5_chr5,
                             wh = 1,
                             geno = chr5_geno,
                             annot = chr5_map,
                             dim = 2)

# ## SV region 2 size 
# North_chr1_map %>% 
#   slice(501:nrow(North_chr1_map)) %>% 
#   summarise(last_pos = last(position), 
#             first_pos = first(position)) %>% 
#   mutate(SV_size = last_pos-first_pos, 
#          SV_size_MB = (last_pos-first_pos)/1000000)
## inversion plots per SNP cluster on the CHR
plot(chr5_inver1)
plot(chr5_inver2)
plot(chr5_inver3) ## SV FOUND
plot(chr5_inver4)
plot(chr5_inver5)


## chr5 sv size 

chr5_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position),
            end = last(position)) %>% 
  mutate(sv_size = end-start, 
         sv_size_mb = sv_size/1000000)

chr5_inver3$datin$y %>% 
  as_tibble() %>% 
  write_csv('chr5_inver3_100SNPS_MDS.csv')

invGenotypes(chr5_inver3) %>% 
  as_tibble() %>% 
  write_csv('chr5_inver3_100SNPS_Inversion_genos.csv')


# plot chr5  sv --------------------------------------------------------------

chr5_mds = read_csv('chr5_inver3_100SNPS_MDS.csv')
chr5_inver_geno = read_csv('chr5_inver3_100SNPS_Inversion_genos.csv')
chr5_ped = read_table2('AC_New_CHRSET_5.ped', 
                       col_names = F) %>% 
  dplyr::select(1:2)

chr5_inver = bind_cols(chr5_ped, 
                       chr5_mds, 
                       chr5_inver_geno)


shades_of_BigD = c('#0583F2', 
                            '#F28705',
                            '#F20530')
                            
# View(North_chr1_Big_D_Energy)
chr5_Big_Plot_Energy = ggplot(data = chr5_inver, 
                                    aes(x = V1, 
                                        y = V2, 
                                        col = value))+
  geom_point()+
  scale_colour_manual(values = shades_of_BigD)+
  # facet_grid(.~label)+
  labs(x = 'MDS axis 1', 
       y = 'MDS axis 2', 
       col = 'Inversion genotypes', 
       title = 'Chr5 11.2 Mb')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

chr5_Big_Plot_Energy

ggsave('Chr5_putate_sv.tiff', 
       plot = chr5_Big_Plot_Energy, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)
# Chr7 --------------------------------------------------------------------
chr7 = read.plink(bed = 'AC_New_CHRSET_7.bed', 
                  bim = 'AC_New_CHRSET_7.bim', 
                  fam = 'AC_New_CHRSET_7.fam')

chr7_geno = chr7$genotypes
chr7_map = chr7$map 
identical(chr7_map[,2], 
          colnames(chr7_geno))


dim(chr7_map)

## SNP cluster 1
chr7_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr7 = data.frame(chr = 7, 
                        LBP = 2047688, 
                        RBP = 24571292, 
                        reg = 'inver1')

chr7_inver1 = invClust(roi = ROI_1_chr7, 
                       wh = 1, 
                       geno = chr7_geno, 
                       annot = chr7_map, 
                       dim = 2)

## SNP cluster 2
chr7_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr7 = data.frame(chr = 7, 
                        LBP = 24571968, 
                        RBP = 32901931, 
                        reg = 'inver2')

chr7_inver2 = invClust(roi = ROI_2_chr7, 
                       wh = 1, 
                       geno = chr7_geno, 
                       annot = chr7_map, 
                       dim = 2)

## SNP cluster 3
chr7_map %>% 
  slice(201:nrow(chr7_map)) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr7 = data.frame(chr = 7, 
                        LBP = 33159264, 
                        RBP = 64728905, 
                        reg = 'inver3')

chr7_inver3 = invClust(roi = ROI_3_chr7, 
                       wh = 1, 
                       geno = chr7_geno, 
                       annot = chr7_map, 
                       dim = 2)


## inversion plots per SNP cluster on the CHR
plot(chr7_inver1)
plot(chr7_inver2) ## sv found
plot(chr7_inver3) 


## chr7 sv size 

chr7_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position),
            end = last(position)) %>% 
  mutate(sv_size = end-start, 
         sv_size_mb = sv_size/1000000)

chr7_inver2$datin$y %>% 
  as_tibble() %>% 
  write_csv('chr7_inver2_100SNPS_MDS.csv')

invGenotypes(chr7_inver2) %>% 
  as_tibble() %>% 
  write_csv('chr7_inver2_100SNPS_Inversion_genos.csv')


# plot chr7  sv --------------------------------------------------------------

chr7_mds = read_csv('chr7_inver2_100SNPS_MDS.csv')
chr7_inver_geno = read_csv('chr7_inver2_100SNPS_Inversion_genos.csv')
chr7_ped = read_table2('AC_New_CHRSET_7.ped', 
                       col_names = F) %>% 
  dplyr::select(1:2)

chr7_inver = bind_cols(chr7_ped, 
                       chr7_mds, 
                       chr7_inver_geno)


shades_of_BigD = c('#0583F2', 
                            '#F28705',
                            '#F20530')
                            
# View(North_chr1_Big_D_Energy)
chr7_Big_Plot_Energy = ggplot(data = chr7_inver, 
                              aes(x = V1, 
                                  y = V2, 
                                  col = value))+
  geom_point()+
  scale_colour_manual(values = shades_of_BigD)+
  # facet_grid(.~label)+
  labs(x = 'MDS axis 1', 
       y = 'MDS axis 2', 
       col = 'Inversion genotypes', 
       title = 'Chr7 15.7 Mb')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

chr7_Big_Plot_Energy

ggsave('chr7_putate_sv.tiff', 
       plot = chr7_Big_Plot_Energy, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)

# Chr9 --------------------------------------------------------------------
chr9 = read.plink(bed = 'AC_New_CHRSET_9.bed', 
                  bim = 'AC_New_CHRSET_9.bim', 
                  fam = 'AC_New_CHRSET_9.fam')

chr9_geno = chr9$genotypes
chr9_map = chr9$map 
identical(chr9_map[,2], 
          colnames(chr9_geno))

dim(chr9_map)

## SNP cluster 1
chr9_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr9 = data.frame(chr = 9, 
                        LBP = 3571127, 
                        RBP = 26036010, 
                        reg = 'inver1')

chr9_inver1 = invClust(roi = ROI_1_chr9, 
                       wh = 1, 
                       geno = chr9_geno, 
                       annot = chr9_map, 
                       dim = 2)

## SNP cluster 2
chr9_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr9 = data.frame(chr = 9, 
                        LBP = 26036012, 
                        RBP = 38401050, 
                        reg = 'inver2')

chr9_inver2 = invClust(roi = ROI_2_chr9, 
                       wh = 1, 
                       geno = chr9_geno, 
                       annot = chr9_map, 
                       dim = 2)

## SNP cluster 3
chr9_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr9 = data.frame(chr = 9, 
                        LBP = 38481505, 
                        RBP = 53111639, 
                        reg = 'inver3')

chr9_inver3 = invClust(roi = ROI_3_chr9, 
                       wh = 1, 
                       geno = chr9_geno, 
                       annot = chr9_map, 
                       dim = 2)

## SNP cluster 4
chr9_map %>% 
  slice(301:400) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_4_chr9 = data.frame(chr = 9,
                        LBP = 53142371,
                        RBP = 64888731,
                        reg = 'inver4')

chr9_inver4 = invClust(roi = ROI_4_chr9,
                       wh = 1,
                       geno = chr9_geno,
                       annot = chr9_map,
                       dim = 2)


## SNP cluster 5
chr9_map %>% 
  slice(401:nrow(chr9_map)) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_5_chr9 = data.frame(chr = 9,
                        LBP = 65323041,
                        RBP = 87135312,
                        reg = 'inver5')

chr9_inver5 = invClust(roi = ROI_5_chr9,
                       wh = 1,
                       geno = chr9_geno,
                       annot = chr9_map,
                       dim = 2)

## inversion plots per SNP cluster on the CHR
plot(chr9_inver1)
plot(chr9_inver2)
plot(chr9_inver3)
plot(chr9_inver4) 
plot(chr9_inver5)


## chr9 sv size 

chr9_map %>% 
  slice(301:400) %>% 
  summarize(start = first(position),
            end = last(position)) %>% 
  mutate(sv_size = end-start, 
         sv_size_mb = sv_size/1000000)

chr9_inver4$datin$y %>% 
  as_tibble() %>% 
  write_csv('chr9_inver4_100SNPS_MDS.csv')

invGenotypes(chr9_inver4) %>% 
  as_tibble() %>% 
  write_csv('chr9_inver4_100SNPS_Inversion_genos.csv')


# plot chr9  sv --------------------------------------------------------------

chr9_mds = read_csv('chr9_inver4_100SNPS_MDS.csv')
chr9_inver_geno = read_csv('chr9_inver4_100SNPS_Inversion_genos.csv')
chr9_ped = read_table2('AC_New_CHRSET_5.ped', 
                       col_names = F) %>% 
  dplyr::select(1:2)

chr9_inver = bind_cols(chr9_ped, 
                       chr9_mds, 
                       chr9_inver_geno)


shades_of_BigD = c('#0583F2', 
                            '#F28705',
                            '#F20530')
                            
# View(North_chr1_Big_D_Energy)
chr9_Big_Plot_Energy = ggplot(data = chr9_inver, 
                              aes(x = V1, 
                                  y = V2, 
                                  col = value))+
  geom_point()+
  scale_colour_manual(values = shades_of_BigD)+
  # facet_grid(.~label)+
  labs(x = 'MDS axis 1', 
       y = 'MDS axis 2', 
       col = 'Inversion genotypes', 
       title = 'chr9 11.2 Mb')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

chr9_Big_Plot_Energy

ggsave('chr9_putate_sv.tiff', 
       plot = chr9_Big_Plot_Energy, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)





# chr 11 ------------------------------------------------------------------
chr11 = read.plink(bed = 'AC_New_CHRSET_11.bed', 
                  bim = 'AC_New_CHRSET_11.bim', 
                  fam = 'AC_New_CHRSET_11.fam')

chr11_geno = chr11$genotypes
chr11_map = chr11$map 
identical(chr11_map[,2], 
          colnames(chr11_geno))

dim(chr11_map)

## SNP cluster 1
chr11_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr11 = data.frame(chr = 11, 
                        LBP = 248694, 
                        RBP = 13240159, 
                        reg = 'inver1')

chr11_inver1 = invClust(roi = ROI_1_chr11, 
                       wh = 1, 
                       geno = chr11_geno, 
                       annot = chr11_map, 
                       dim = 2)

## SNP cluster 2
chr11_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr11 = data.frame(chr = 11, 
                        LBP = 13260402, 
                        RBP = 26892078, 
                        reg = 'inver2')

chr11_inver2 = invClust(roi = ROI_2_chr11, 
                       wh = 1, 
                       geno = chr11_geno, 
                       annot = chr11_map, 
                       dim = 2)

## SNP cluster 3
chr11_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr11 = data.frame(chr = 11, 
                        LBP = 26946851, 
                        RBP = 49959767, 
                        reg = 'inver3')

chr11_inver3 = invClust(roi = ROI_3_chr11, 
                       wh = 1, 
                       geno = chr11_geno, 
                       annot = chr11_map, 
                       dim = 2)

## SNP cluster 4
chr11_map %>% 
  slice(301:400) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_4_chr11 = data.frame(chr = 11,
                        LBP = 49960044,
                        RBP = 53713982,
                        reg = 'inver4')

chr11_inver4 = invClust(roi = ROI_4_chr11,
                       wh = 1,
                       geno = chr11_geno,
                       annot = chr11_map,
                       dim = 2)


## SNP cluster 5
chr11_map %>% 
  slice(401:501) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_5_chr11 = data.frame(chr = 11,
                        LBP = 53716090,
                        RBP = 59171961,
                        reg = 'inver5')

chr11_inver5 = invClust(roi = ROI_5_chr11,
                       wh = 1,
                       geno = chr11_geno,
                       annot = chr11_map,
                       dim = 2)

## SNP cluster 6
chr11_map %>% 
  slice(501:nrow(chr11_map)) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_6_chr11 = data.frame(chr = 11,
                         LBP = 59171961,
                         RBP = 87346350,
                         reg = 'inver5')

chr11_inver6 = invClust(roi = ROI_6_chr11,
                        wh = 1,
                        geno = chr11_geno,
                        annot = chr11_map,
                        dim = 2)

## inversion plots per SNP cluster on the CHR
plot(chr11_inver1)
plot(chr11_inver2) 
plot(chr11_inver3) ## Only inversions
plot(chr11_inver4) 
plot(chr11_inver5) 
plot(chr11_inver6)


## chr11 sv size 

# chr11_map %>% 
#   slice(101:200) %>% 
#   summarize(start = first(position),
#             end = last(position)) %>% 
#   mutate(sv_size = end-start, 
#          sv_size_mb = sv_size/1000000)
# 
# chr11_inver2$datin$y %>% 
#   as_tibble() %>% 
#   write_csv('chr11_inver2_100SNPS_MDS.csv')
# 
# invGenotypes(chr11_inver2) %>% 
#   as_tibble() %>% 
#   write_csv('chr11_inver2_100SNPS_Inversion_genos.csv')

chr11_map %>%
  slice(201:300) %>%
  summarize(start = first(position),
            end = last(position)) %>%
  mutate(sv_size = end-start,
         sv_size_mb = sv_size/1000000)

chr11_inver3$datin$y %>%
  as_tibble() %>%
  write_csv('chr11_inver3_100SNPS_MDS.csv')

invGenotypes(chr11_inver3) %>%
  as_tibble() %>%
  write_csv('chr11_inver3_100SNPS_Inversion_genos.csv')
# 
# 
# chr11_map %>% 
#   slice(301:400) %>% 
#   summarize(start = first(position),
#             end = last(position)) %>% 
#   mutate(sv_size = end-start, 
#          sv_size_mb = sv_size/1000000)
# 
# chr11_inver4$datin$y %>% 
#   as_tibble() %>% 
#   write_csv('chr11_inver4_100SNPS_MDS.csv')
# 
# invGenotypes(chr11_inver4) %>% 
#   as_tibble() %>% 
#   write_csv('chr11_inver4_100SNPS_Inversion_genos.csv')
# 
# 
# chr11_map %>% 
#   slice(401:500) %>% 
#   summarize(start = first(position),
#             end = last(position)) %>% 
#   mutate(sv_size = end-start, 
#          sv_size_mb = sv_size/1000000)
# 
# chr11_inver5$datin$y %>% 
#   as_tibble() %>% 
#   write_csv('chr11_inver5_100SNPS_MDS.csv')
# 
# invGenotypes(chr11_inver5) %>% 
#   as_tibble() %>% 
#   write_csv('chr11_inver5_100SNPS_Inversion_genos.csv')
# 

# chr11 data formatting ---------------------------------------------------


chr11_ped = read_table2('AC_New_CHRSET_5.ped', 
                        col_names = F) %>% 
  dplyr::select(1:2)

# chr11_mds2 = read_csv('chr11_inver2_100SNPS_MDS.csv')
# chr11_inver_geno2 = read_csv('chr11_inver2_100SNPS_Inversion_genos.csv')
# 
# chr11_inver2 = bind_cols(chr11_ped, 
#                        chr11_mds2, 
#                        chr11_inver_geno2)

# label = rep('inversion1', 
#             nrow(chr11_inver2))
# chr11_inver2 = bind_cols(chr11_inver2, 
#                          label)
# 
# 
chr11_mds3 = read_csv('chr11_inver3_100SNPS_MDS.csv')
chr11_inver_geno3 = read_csv('chr11_inver3_100SNPS_Inversion_genos.csv')

chr11_inver3 = bind_cols(chr11_ped,
                         chr11_mds3,
                         chr11_inver_geno3)
# 
# label = rep('inversion2', 
#             nrow(chr11_inver3))
# chr11_inver3 = bind_cols(chr11_inver3, 
#                          label)
# 
# 
# chr11_mds4 = read_csv('chr11_inver4_100SNPS_MDS.csv')
# chr11_inver_geno4 = read_csv('chr11_inver4_100SNPS_Inversion_genos.csv')
# 
# chr11_inver4 = bind_cols(chr11_ped, 
#                          chr11_mds4, 
#                          chr11_inver_geno4)
# 
# label = rep('inversion3', 
#             nrow(chr11_inver4))
# chr11_inver4 = bind_cols(chr11_inver4, 
#                          label)
# 
# 
# 
# chr11_mds5 = read_csv('chr11_inver5_100SNPS_MDS.csv')
# chr11_inver_geno5 = read_csv('chr11_inver5_100SNPS_Inversion_genos.csv')
# 
# chr11_inver5 = bind_cols(chr11_ped, 
#                          chr11_mds2, 
#                          chr11_inver_geno2)
# 
# label = rep('inversion4', 
#             nrow(chr11_inver5))
# chr11_inver5 = bind_cols(chr11_inver5, 
#                          label)
# 
# 
# chr11_all_sv = bind_rows(chr11_inver2, 
#                          chr11_inver3, 
#                          chr11_inver4, 
#                          chr11_inver5) %>% 
#   rename(inversion = ...9)

# plot chr11  sv --------------------------------------------------------------

shades_of_BigD = c('#0583F2', 
                            '#F28705',
                            '#F20530')
                            
# View(North_chr1_Big_D_Energy)
chr11_Big_Plot_Energy = ggplot(data = chr11_inver3, 
                              aes(x = V1, 
                                  y = V2, 
                                  col = value))+
  geom_point()+
  scale_colour_manual(values = shades_of_BigD)+
  # facet_grid(.~inversion)+
  labs(x = 'MDS axis 1', 
       y = 'MDS axis 2', 
       col = 'Inversion genotypes', 
       title = 'chr11 23.0 Mb')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

chr11_Big_Plot_Energy

ggsave('chr11_putate_sv.tiff', 
       plot = chr11_Big_Plot_Energy, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)








# chr 12 ------------------------------------------------------------------

chr12 = read.plink(bed = 'AC_New_CHRSET_12.bed', 
                  bim = 'AC_New_CHRSET_12.bim', 
                  fam = 'AC_New_CHRSET_12.fam')

chr12_geno = chr12$genotypes
chr12_map = chr12$map 
identical(chr12_map[,2], 
          colnames(chr12_geno))

dim(chr12_map)

## SNP cluster 1
chr12_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr12 = data.frame(chr = 12, 
                        LBP = 4943808, 
                        RBP = 16976261, 
                        reg = 'inver1')

chr12_inver1 = invClust(roi = ROI_1_chr12, 
                       wh = 1, 
                       geno = chr12_geno, 
                       annot = chr12_map, 
                       dim = 2)

## SNP cluster 2
chr12_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr12 = data.frame(chr = 12, 
                        LBP = 17149318, 
                        RBP = 32764904, 
                        reg = 'inver2')

chr12_inver2 = invClust(roi = ROI_2_chr12, 
                       wh = 1, 
                       geno = chr12_geno, 
                       annot = chr12_map, 
                       dim = 2)

## SNP cluster 3
chr12_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr12 = data.frame(chr = 12, 
                        LBP = 32966855, 
                        RBP = 46370721, 
                        reg = 'inver3')

chr12_inver3 = invClust(roi = ROI_3_chr12, 
                       wh = 1, 
                       geno = chr12_geno, 
                       annot = chr12_map, 
                       dim = 2)

## SNP cluster 4
chr12_map %>% 
  slice(301:nrow(chr12_map)) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_4_chr12 = data.frame(chr = 11,
                        LBP = 46701483,
                        RBP = 56958464,
                        reg = 'inver4')

chr12_inver4 = invClust(roi = ROI_4_chr12,
                       wh = 1,
                       geno = chr12_geno,
                       annot = chr12_map,
                       dim = 2)


## inversion plots per SNP cluster on the CHR
plot(chr12_inver1)
plot(chr12_inver2)
plot(chr12_inver3)
plot(chr12_inver4) ## maybe?


## chr12 sv size 

chr12_map %>% 
  slice(301:nrow(chr12_map)) %>% 
  summarize(start = first(position),
            end = last(position)) %>% 
  mutate(sv_size = end-start, 
         sv_size_mb = sv_size/1000000)

chr12_inver4$datin$y %>% 
  as_tibble() %>% 
  write_csv('chr12_inver4_100SNPS_MDS.csv')

invGenotypes(chr12_inver4) %>% 
  as_tibble() %>% 
  write_csv('chr12_inver4_100SNPS_Inversion_genos.csv')


# plot chr12  sv --------------------------------------------------------------

chr12_mds = read_csv('chr12_inver4_100SNPS_MDS.csv')
chr12_inver_geno = read_csv('chr12_inver4_100SNPS_Inversion_genos.csv')
chr12_ped = read_table2('AC_New_CHRSET_5.ped', 
                       col_names = F) %>% 
  dplyr::select(1:2)

chr12_inver = bind_cols(chr12_ped, 
                       chr12_mds, 
                       chr12_inver_geno)


shades_of_BigD = c('#0583F2', 
                            '#F28705',
                            '#F20530')
                            
# View(North_chr1_Big_D_Energy)
chr12_Big_Plot_Energy = ggplot(data = chr12_inver, 
                              aes(x = V1, 
                                  y = V2, 
                                  col = value))+
  geom_point()+
  scale_colour_manual(values = shades_of_BigD)+
  # facet_grid(.~label)+
  labs(x = 'MDS axis 1', 
       y = 'MDS axis 2', 
       col = 'Inversion genotypes', 
       title = 'chr12 10.3 Mb')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

chr12_Big_Plot_Energy

ggsave('chr12_putate_sv.tiff', 
       plot = chr12_Big_Plot_Energy, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)





chr22 = read.plink(bed = 'AC_New_CHRSET_22.bed', 
                  bim = 'AC_New_CHRSET_22.bim', 
                  fam = 'AC_New_CHRSET_22.fam')

chr22_geno = chr22$genotypes
chr22_map = chr22$map 
identical(chr22_map[,2], 
          colnames(chr22_geno))

dim(chr22_map)

## SNP cluster 1
chr22_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr22 = data.frame(chr = 22, 
                        LBP = 112929, 
                        RBP = 11498960, 
                        reg = 'inver1')

chr22_inver1 = invClust(roi = ROI_1_chr22, 
                       wh = 1, 
                       geno = chr22_geno, 
                       annot = chr22_map, 
                       dim = 2)

## SNP cluster 2
chr22_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr22 = data.frame(chr = 22, 
                        LBP = 11575805, 
                        RBP = 30886025, 
                        reg = 'inver2')

chr22_inver2 = invClust(roi = ROI_2_chr22, 
                       wh = 1, 
                       geno = chr22_geno, 
                       annot = chr22_map, 
                       dim = 2)

## SNP cluster 3
chr22_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr22 = data.frame(chr = 22, 
                        LBP = 31347076, 
                        RBP = 44816796, 
                        reg = 'inver3')

chr22_inver3 = invClust(roi = ROI_3_chr22, 
                       wh = 1, 
                       geno = chr22_geno, 
                       annot = chr22_map, 
                       dim = 2)

## SNP cluster 4
chr22_map %>% 
  slice(301:nrow(chr22_map)) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_4_chr22 = data.frame(chr = 22,
                        LBP = 45121392,
                        RBP = 54338362,
                        reg = 'inver4')

chr22_inver4 = invClust(roi = ROI_4_chr22,
                       wh = 1,
                       geno = chr22_geno,
                       annot = chr22_map,
                       dim = 2)



## inversion plots per SNP cluster on the CHR
plot(chr22_inver1)
plot(chr22_inver2)
plot(chr22_inver3)
plot(chr22_inver4) ## maybe


## chr22 sv size 

chr22_map %>% 
  slice(301:nrow(chr22_map)) %>% 
  summarize(start = first(position),
            end = last(position)) %>% 
  mutate(sv_size = end-start, 
         sv_size_mb = sv_size/1000000)

chr22_inver4$datin$y %>% 
  as_tibble() %>% 
  write_csv('chr22_inver4_100SNPS_MDS.csv')

invGenotypes(chr22_inver4) %>% 
  as_tibble() %>% 
  write_csv('chr22_inver4_100SNPS_Inversion_genos.csv')


# plot chr22  sv --------------------------------------------------------------

chr22_mds = read_csv('chr22_inver4_100SNPS_MDS.csv')
chr22_inver_geno = read_csv('chr22_inver4_100SNPS_Inversion_genos.csv')
chr22_ped = read_table2('AC_New_CHRSET_5.ped', 
                       col_names = F) %>% 
  dplyr::select(1:2)

chr22_inver = bind_cols(chr22_ped, 
                       chr22_mds, 
                       chr22_inver_geno)


shades_of_BigD = c('#0583F2', 
                            '#F28705',
                            '#F20530')
                            
# View(North_chr1_Big_D_Energy)
chr22_Big_Plot_Energy = ggplot(data = chr22_inver, 
                              aes(x = V1, 
                                  y = V2, 
                                  col = value))+
  geom_point()+
  scale_colour_manual(values = shades_of_BigD)+
  # facet_grid(.~label)+
  labs(x = 'MDS axis 1', 
       y = 'MDS axis 2', 
       col = 'Inversion genotypes', 
       title = 'chr22 9.2 Mb')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

chr22_Big_Plot_Energy

ggsave('chr22_putate_sv.tiff', 
       plot = chr22_Big_Plot_Energy, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)




# chr 23 ------------------------------------------------------------------


chr23 = read.plink(bed = 'AC_New_CHRSET_23.bed', 
                  bim = 'AC_New_CHRSET_23.bim', 
                  fam = 'AC_New_CHRSET_23.fam')

chr23_geno = chr23$genotypes
chr23_map = chr23$map 
identical(chr23_map[,2], 
          colnames(chr23_geno))

dim(chr23_map)

## SNP cluster 1
chr23_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr23 = data.frame(chr = 23, 
                        LBP = 851590, 
                        RBP = 10742161, 
                        reg = 'inver1')

chr23_inver1 = invClust(roi = ROI_1_chr23, 
                       wh = 1, 
                       geno = chr23_geno, 
                       annot = chr23_map, 
                       dim = 2)

## SNP cluster 2
chr23_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr23 = data.frame(chr = 23, 
                        LBP = 10821290, 
                        RBP = 22240440, 
                        reg = 'inver2')

chr23_inver2 = invClust(roi = ROI_2_chr23, 
                       wh = 1, 
                       geno = chr23_geno, 
                       annot = chr23_map, 
                       dim = 2)

## SNP cluster 3
chr23_map %>% 
  slice(201:nrow(chr23_map)) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr23 = data.frame(chr = 23, 
                        LBP = 22254968, 
                        RBP = 31713869, 
                        reg = 'inver3')

chr23_inver3 = invClust(roi = ROI_3_chr23, 
                       wh = 1, 
                       geno = chr23_geno, 
                       annot = chr23_map, 
                       dim = 2)


## inversion plots per SNP cluster on the CHR
plot(chr23_inver1)
plot(chr23_inver2)
plot(chr23_inver3)## maybe


## chr23 sv size 

chr23_map %>% 
  slice(201:nrow(chr23_map)) %>% 
  summarize(start = first(position),
            end = last(position)) %>% 
  mutate(sv_size = end-start, 
         sv_size_mb = sv_size/1000000)

chr23_inver3$datin$y %>% 
  as_tibble() %>% 
  write_csv('chr23_inver3_100SNPS_MDS.csv')

invGenotypes(chr23_inver3) %>% 
  as_tibble() %>% 
  write_csv('chr23_inver3_100SNPS_Inversion_genos.csv')


# plot chr23  sv --------------------------------------------------------------

chr23_mds = read_csv('chr23_inver3_100SNPS_MDS.csv')
chr23_inver_geno = read_csv('chr23_inver3_100SNPS_Inversion_genos.csv')
chr23_ped = read_table2('AC_New_CHRSET_5.ped', 
                       col_names = F) %>% 
  dplyr::select(1:2)

chr23_inver = bind_cols(chr23_ped, 
                       chr23_mds, 
                       chr23_inver_geno)


shades_of_BigD = c('#0583F2', 
                            '#F28705',
                            '#F20530')
                            
# View(North_chr1_Big_D_Energy)
chr23_Big_Plot_Energy = ggplot(data = chr23_inver, 
                              aes(x = V1, 
                                  y = V2, 
                                  col = value))+
  geom_point()+
  scale_colour_manual(values = shades_of_BigD)+
  # facet_grid(.~label)+
  labs(x = 'MDS axis 1', 
       y = 'MDS axis 2', 
       col = 'Inversion genotypes', 
       title = 'chr23 9.5 Mb')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

chr23_Big_Plot_Energy

ggsave('chr23_putate_sv.tiff', 
       plot = chr23_Big_Plot_Energy, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)






# chr 26 ------------------------------------------------------------------
chr26 = read.plink(bed = 'AC_New_CHRSET_26.bed', 
                  bim = 'AC_New_CHRSET_26.bim', 
                  fam = 'AC_New_CHRSET_26.fam')

chr26_geno = chr26$genotypes
chr26_map = chr26$map 
identical(chr26_map[,2], 
          colnames(chr26_geno))

dim(chr26_map)

## SNP cluster 1
chr26_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr26 = data.frame(chr = 26, 
                        LBP = 883499, 
                        RBP = 11053831, 
                        reg = 'inver1')

chr26_inver1 = invClust(roi = ROI_1_chr26, 
                       wh = 1, 
                       geno = chr26_geno, 
                       annot = chr26_map, 
                       dim = 2)

## SNP cluster 2
chr26_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr26 = data.frame(chr = 26, 
                        LBP = 11221098, 
                        RBP = 25552361, 
                        reg = 'inver2')

chr26_inver2 = invClust(roi = ROI_2_chr26, 
                       wh = 1, 
                       geno = chr26_geno, 
                       annot = chr26_map, 
                       dim = 2)

## SNP cluster 3
chr26_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr26 = data.frame(chr = 26, 
                        LBP = 25619512, 
                        RBP = 29581331, 
                        reg = 'inver3')

chr26_inver3 = invClust(roi = ROI_3_chr26, 
                       wh = 1, 
                       geno = chr26_geno, 
                       annot = chr26_map, 
                       dim = 2)

## SNP cluster 4
chr26_map %>% 
  slice(301:nrow(chr26_map)) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_4_chr26 = data.frame(chr = 26,
                        LBP = 29599392,
                        RBP = 39785897,
                        reg = 'inver4')

chr26_inver4 = invClust(roi = ROI_4_chr26,
                       wh = 1,
                       geno = chr26_geno,
                       annot = chr26_map,
                       dim = 2)



## inversion plots per SNP cluster on the CHR
plot(chr26_inver1)
plot(chr26_inver2)
plot(chr26_inver3)
plot(chr26_inver4) 


## chr26 sv size 

chr26_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position),
            end = last(position)) %>% 
  mutate(sv_size = end-start, 
         sv_size_mb = sv_size/1000000)

chr26_inver3$datin$y %>% 
  as_tibble() %>% 
  write_csv('chr26_inver3_100SNPS_MDS.csv')

invGenotypes(chr26_inver3) %>% 
  as_tibble() %>% 
  write_csv('chr26_inver3_100SNPS_Inversion_genos.csv')


# plot chr26  sv --------------------------------------------------------------

chr26_mds = read_csv('chr26_inver3_100SNPS_MDS.csv')
chr26_inver_geno = read_csv('chr26_inver3_100SNPS_Inversion_genos.csv')
chr26_ped = read_table2('AC_New_CHRSET_5.ped', 
                       col_names = F) %>% 
  dplyr::select(1:2)

chr26_inver = bind_cols(chr26_ped, 
                       chr26_mds, 
                       chr26_inver_geno)


shades_of_BigD = c('#0583F2', 
                            '#F28705',
                            '#F20530')
                            
# View(North_chr1_Big_D_Energy)
chr26_Big_Plot_Energy = ggplot(data = chr26_inver, 
                              aes(x = V1, 
                                  y = V2, 
                                  col = value))+
  geom_point()+
  scale_colour_manual(values = shades_of_BigD)+
  # facet_grid(.~label)+
  labs(x = 'MDS axis 1', 
       y = 'MDS axis 2', 
       col = 'Inversion genotypes', 
       title = 'chr26 11.2 Mb')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

chr26_Big_Plot_Energy

ggsave('chr26_putate_sv.tiff', 
       plot = chr26_Big_Plot_Energy, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)

# chr27 -------------------------------------------------------------------

chr27 = read.plink(bed = 'AC_New_CHRSET_27.bed', 
                  bim = 'AC_New_CHRSET_27.bim', 
                  fam = 'AC_New_CHRSET_27.fam')

chr27_geno = chr27$genotypes
chr27_map = chr27$map 
identical(chr27_map[,2], 
          colnames(chr27_geno))

dim(chr27_map)

## SNP cluster 1
chr27_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr27 = data.frame(chr = 27, 
                        LBP = 7427836, 
                        RBP = 26082025, 
                        reg = 'inver1')

chr27_inver1 = invClust(roi = ROI_1_chr27, 
                       wh = 1, 
                       geno = chr27_geno, 
                       annot = chr27_map, 
                       dim = 2)

## SNP cluster 2
chr27_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr27 = data.frame(chr = 27, 
                        LBP = 26256039, 
                        RBP = 36620257, 
                        reg = 'inver2')

chr27_inver2 = invClust(roi = ROI_2_chr27, 
                       wh = 1, 
                       geno = chr27_geno, 
                       annot = chr27_map, 
                       dim = 2)

## SNP cluster 3
chr27_map %>% 
  slice(201:nrow(chr27_map)) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr27 = data.frame(chr = 27, 
                        LBP = 36706140, 
                        RBP = 45228603, 
                        reg = 'inver3')

chr27_inver3 = invClust(roi = ROI_3_chr27, 
                       wh = 1, 
                       geno = chr27_geno, 
                       annot = chr27_map, 
                       dim = 2)

## inversion plots per SNP cluster on the CHR
plot(chr27_inver1) #yep
plot(chr27_inver2)
plot(chr27_inver3)


## chr27 sv size 

chr27_map %>% 
  slice(1:100) %>% 
  summarize(start = first(position),
            end = last(position)) %>% 
  mutate(sv_size = end-start, 
         sv_size_mb = sv_size/1000000)

chr27_inver1$datin$y %>% 
  as_tibble() %>% 
  write_csv('chr27_inver1_100SNPS_MDS.csv')

invGenotypes(chr27_inver1) %>% 
  as_tibble() %>% 
  write_csv('chr27_inver1_100SNPS_Inversion_genos.csv')


# plot chr27  sv --------------------------------------------------------------

chr27_mds = read_csv('chr27_inver1_100SNPS_MDS.csv')
chr27_inver_geno = read_csv('chr27_inver1_100SNPS_Inversion_genos.csv')
chr27_ped = read_table2('AC_New_CHRSET_5.ped', 
                       col_names = F) %>% 
  dplyr::select(1:2)

chr27_inver = bind_cols(chr27_ped, 
                       chr27_mds, 
                       chr27_inver_geno)


shades_of_BigD = c('#0583F2', 
                            '#F28705',
                            '#F20530')
                            
# View(North_chr1_Big_D_Energy)
chr27_Big_Plot_Energy = ggplot(data = chr27_inver, 
                              aes(x = V1, 
                                  y = V2, 
                                  col = value))+
  geom_point()+
  scale_colour_manual(values = shades_of_BigD)+
  # facet_grid(.~label)+
  labs(x = 'MDS axis 1', 
       y = 'MDS axis 2', 
       col = 'Inversion genotypes', 
       title = 'chr27 11.2 Mb')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

chr27_Big_Plot_Energy

ggsave('chr27_putate_sv.tiff', 
       plot = chr27_Big_Plot_Energy, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)

# chr30 -------------------------------------------------------------------

chr30 = read.plink(bed = 'AC_New_CHRSET_30.bed', 
                  bim = 'AC_New_CHRSET_30.bim', 
                  fam = 'AC_New_CHRSET_30.fam')

chr30_geno = chr30$genotypes
chr30_map = chr30$map 
identical(chr30_map[,2], 
          colnames(chr30_geno))

dim(chr30_map)

## SNP cluster 1
chr30_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr30 = data.frame(chr = 30, 
                        LBP = 867959, 
                        RBP = 23790862, 
                        reg = 'inver1')

chr30_inver1 = invClust(roi = ROI_1_chr30, 
                       wh = 1, 
                       geno = chr30_geno, 
                       annot = chr30_map, 
                       dim = 2)

## SNP cluster 2
chr30_map %>% 
  slice(101:200) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr30 = data.frame(chr = 30, 
                        LBP = 23803357, 
                        RBP = 32694308, 
                        reg = 'inver2')

chr30_inver2 = invClust(roi = ROI_2_chr30, 
                       wh = 1, 
                       geno = chr30_geno, 
                       annot = chr30_map, 
                       dim = 2)

## SNP cluster 3
chr30_map %>% 
  slice(201:300) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_3_chr30 = data.frame(chr = 30, 
                        LBP = 32694374, 
                        RBP = 43436389, 
                        reg = 'inver3')

chr30_inver3 = invClust(roi = ROI_3_chr30, 
                       wh = 1, 
                       geno = chr30_geno, 
                       annot = chr30_map, 
                       dim = 2)

## SNP cluster 4
chr30_map %>% 
  slice(301:nrow(chr30_map)) %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_4_chr30 = data.frame(chr = 30,
                        LBP = 43513470,
                        RBP = 50061167,
                        reg = 'inver4')

chr30_inver4 = invClust(roi = ROI_4_chr30,
                       wh = 1,
                       geno = chr30_geno,
                       annot = chr30_map,
                       dim = 2)



## inversion plots per SNP cluster on the CHR
plot(chr30_inver1)
plot(chr30_inver2)
plot(chr30_inver3)
plot(chr30_inver4)


# ## chr30 sv size 
# 
# chr30_map %>% 
#   slice(301:400) %>% 
#   summarize(start = first(position),
#             end = last(position)) %>% 
#   mutate(sv_size = end-start, 
#          sv_size_mb = sv_size/1000000)
# 
# chr30_inver4$datin$y %>% 
#   as_tibble() %>% 
#   write_csv('chr30_inver4_100SNPS_MDS.csv')
# 
# invGenotypes(chr30_inver4) %>% 
#   as_tibble() %>% 
#   write_csv('chr30_inver4_100SNPS_Inversion_genos.csv')
# 
# 
# # plot chr30  sv --------------------------------------------------------------
# 
# chr30_mds = read_csv('chr30_inver4_100SNPS_MDS.csv')
# chr30_inver_geno = read_csv('chr30_inver4_100SNPS_Inversion_genos.csv')
# chr30_ped = read_table2('AC_New_CHRSET_5.ped', 
#                        col_names = F) %>% 
#   dplyr::select(1:2)
# 
# chr30_inver = bind_cols(chr30_ped, 
#                        chr30_mds, 
#                        chr30_inver_geno)
# 
# 
# shades_of_BigD = c('#0583F2', 
#                             '#F28705',
#                             '#F20530')
#                             
# # View(North_chr1_Big_D_Energy)
# chr30_Big_Plot_Energy = ggplot(data = chr30_inver, 
#                               aes(x = V1, 
#                                   y = V2, 
#                                   col = value))+
#   geom_point()+
#   scale_colour_manual(values = shades_of_BigD)+
#   # facet_grid(.~label)+
#   labs(x = 'MDS axis 1', 
#        y = 'MDS axis 2', 
#        col = 'Inversion genotypes', 
#        title = 'chr30 11.2 Mb')+
#   theme(panel.grid = element_blank(), 
#         axis.title = element_text(size = 14), 
#         axis.text.y = element_text(size = 12), 
#         axis.text.x = element_text(size = 10, 
#                                    angle = 45, 
#                                    hjust = 1, 
#                                    vjust = 1),
#         legend.title = element_text(size = 14), 
#         legend.text = element_text(size = 12), 
#         strip.background = element_rect(fill = 'white'),
#         strip.text = element_text(size = 12, 
#                                   face = 'bold'))
# 
# chr30_Big_Plot_Energy
# 
# ggsave('chr30_putate_sv.tiff', 
#        plot = chr30_Big_Plot_Energy, 
#        dpi = 'retina', 
#        units = 'cm', 
#        width = 20, 
#        height = 15)
# 


# chr34 -------------------------------------------------------------------

chr34 = read.plink(bed = 'AC_New_CHRSET_34.bed', 
                  bim = 'AC_New_CHRSET_34.bim', 
                  fam = 'AC_New_CHRSET_34.fam')

chr34_geno = chr34$genotypes
chr34_map = chr34$map 
identical(chr34_map[,2], 
          colnames(chr34_geno))

dim(chr34_map)

## SNP cluster 1
chr34_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr34 = data.frame(chr = 34, 
                        LBP = 181141, 
                        RBP = 12002312, 
                        reg = 'inver1')

chr34_inver1 = invClust(roi = ROI_1_chr34, 
                       wh = 1, 
                       geno = chr34_geno, 
                       annot = chr34_map, 
                       dim = 2)

## SNP cluster 2
chr34_map %>% 
  slice(101:nrow(chr34_map)) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr34 = data.frame(chr = 34, 
                        LBP = 12009205, 
                        RBP = 34847981, 
                        reg = 'inver2')

chr34_inver2 = invClust(roi = ROI_2_chr34, 
                       wh = 1, 
                       geno = chr34_geno, 
                       annot = chr34_map, 
                       dim = 2)


## inversion plots per SNP cluster on the CHR
plot(chr34_inver1)
plot(chr34_inver2)


# ## chr34 sv size 
# 
# chr34_map %>% 
#   slice(301:400) %>% 
#   summarize(start = first(position),
#             end = last(position)) %>% 
#   mutate(sv_size = end-start, 
#          sv_size_mb = sv_size/1000000)
# 
# chr34_inver4$datin$y %>% 
#   as_tibble() %>% 
#   write_csv('chr34_inver4_100SNPS_MDS.csv')
# 
# invGenotypes(chr34_inver4) %>% 
#   as_tibble() %>% 
#   write_csv('chr34_inver4_100SNPS_Inversion_genos.csv')
# 
# 
# # plot chr34  sv --------------------------------------------------------------
# 
# chr34_mds = read_csv('chr34_inver4_100SNPS_MDS.csv')
# chr34_inver_geno = read_csv('chr34_inver4_100SNPS_Inversion_genos.csv')
# chr34_ped = read_table2('AC_New_CHRSET_5.ped', 
#                        col_names = F) %>% 
#   dplyr::select(1:2)
# 
# chr34_inver = bind_cols(chr34_ped, 
#                        chr34_mds, 
#                        chr34_inver_geno)
# 
# 
# shades_of_BigD = c('#0583F2', 
#                             '#F28705',
#                             '#F20530')
#                             
# # View(North_chr1_Big_D_Energy)
# chr34_Big_Plot_Energy = ggplot(data = chr34_inver, 
#                               aes(x = V1, 
#                                   y = V2, 
#                                   col = value))+
#   geom_point()+
#   scale_colour_manual(values = shades_of_BigD)+
#   # facet_grid(.~label)+
#   labs(x = 'MDS axis 1', 
#        y = 'MDS axis 2', 
#        col = 'Inversion genotypes', 
#        title = 'chr34 11.2 Mb')+
#   theme(panel.grid = element_blank(), 
#         axis.title = element_text(size = 14), 
#         axis.text.y = element_text(size = 12), 
#         axis.text.x = element_text(size = 10, 
#                                    angle = 45, 
#                                    hjust = 1, 
#                                    vjust = 1),
#         legend.title = element_text(size = 14), 
#         legend.text = element_text(size = 12), 
#         strip.background = element_rect(fill = 'white'),
#         strip.text = element_text(size = 12, 
#                                   face = 'bold'))
# 
# chr34_Big_Plot_Energy
# 
# ggsave('chr34_putate_sv.tiff', 
#        plot = chr34_Big_Plot_Energy, 
#        dpi = 'retina', 
#        units = 'cm', 
#        width = 20, 
#        height = 15)


# chr35 -------------------------------------------------------------------

chr35 = read.plink(bed = 'AC_New_CHRSET_35.bed', 
                  bim = 'AC_New_CHRSET_35.bim', 
                  fam = 'AC_New_CHRSET_35.fam')

chr35_geno = chr35$genotypes
chr35_map = chr35$map 
identical(chr35_map[,2], 
          colnames(chr35_geno))

dim(chr35_map)

## SNP cluster 1
chr35_map %>% 
  slice(1:100) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr35 = data.frame(chr = 35, 
                        LBP = 610665, 
                        RBP = 13407556, 
                        reg = 'inver1')

chr35_inver1 = invClust(roi = ROI_1_chr35, 
                       wh = 1, 
                       geno = chr35_geno, 
                       annot = chr35_map, 
                       dim = 2)

## SNP cluster 2
chr35_map %>% 
  slice(101:nrow(chr35_map)) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_2_chr35 = data.frame(chr = 35, 
                        LBP = 13455118, 
                        RBP = 20652526, 
                        reg = 'inver2')

chr35_inver2 = invClust(roi = ROI_2_chr35, 
                       wh = 1, 
                       geno = chr35_geno, 
                       annot = chr35_map, 
                       dim = 2)


## inversion plots per SNP cluster on the CHR
plot(chr35_inver1)
plot(chr35_inver2) ## maybe

# ## chr35 sv size 
# 
# chr35_map %>% 
#   slice(301:400) %>% 
#   summarize(start = first(position),
#             end = last(position)) %>% 
#   mutate(sv_size = end-start, 
#          sv_size_mb = sv_size/1000000)
# 
# chr35_inver4$datin$y %>% 
#   as_tibble() %>% 
#   write_csv('chr35_inver4_100SNPS_MDS.csv')
# 
# invGenotypes(chr35_inver4) %>% 
#   as_tibble() %>% 
#   write_csv('chr35_inver4_100SNPS_Inversion_genos.csv')
# 
# 
# # plot chr35  sv --------------------------------------------------------------
# 
# chr35_mds = read_csv('chr35_inver4_100SNPS_MDS.csv')
# chr35_inver_geno = read_csv('chr35_inver4_100SNPS_Inversion_genos.csv')
# chr35_ped = read_table2('AC_New_CHRSET_5.ped', 
#                        col_names = F) %>% 
#   dplyr::select(1:2)
# 
# chr35_inver = bind_cols(chr35_ped, 
#                        chr35_mds, 
#                        chr35_inver_geno)
# 
# 
# shades_of_BigD = c('#0583F2', 
#                             '#F28705',
#                             '#F20530')
#                             
# # View(North_chr1_Big_D_Energy)
# chr35_Big_Plot_Energy = ggplot(data = chr35_inver, 
#                               aes(x = V1, 
#                                   y = V2, 
#                                   col = value))+
#   geom_point()+
#   scale_colour_manual(values = shades_of_BigD)+
#   # facet_grid(.~label)+
#   labs(x = 'MDS axis 1', 
#        y = 'MDS axis 2', 
#        col = 'Inversion genotypes', 
#        title = 'chr35 11.2 Mb')+
#   theme(panel.grid = element_blank(), 
#         axis.title = element_text(size = 14), 
#         axis.text.y = element_text(size = 12), 
#         axis.text.x = element_text(size = 10, 
#                                    angle = 45, 
#                                    hjust = 1, 
#                                    vjust = 1),
#         legend.title = element_text(size = 14), 
#         legend.text = element_text(size = 12), 
#         strip.background = element_rect(fill = 'white'),
#         strip.text = element_text(size = 12, 
#                                   face = 'bold'))
# 
# chr35_Big_Plot_Energy
# 
# ggsave('chr35_putate_sv.tiff', 
#        plot = chr35_Big_Plot_Energy, 
#        dpi = 'retina', 
#        units = 'cm', 
#        width = 20, 
#        height = 15)


# chr36 -------------------------------------------------------------------

chr36 = read.plink(bed = 'AC_New_CHRSET_36.bed', 
                  bim = 'AC_New_CHRSET_36.bim', 
                  fam = 'AC_New_CHRSET_36.fam')

chr36_geno = chr36$genotypes
chr36_map = chr36$map 
identical(chr36_map[,2], 
          colnames(chr36_geno))

dim(chr36_map)

## SNP cluster 1
chr36_map %>% 
  slice(1:nrow(chr36_map)) %>% 
  as_tibble() %>% 
  summarize(start = first(position), 
            end = last(position))

ROI_1_chr36 = data.frame(chr = 36, 
                        LBP = 714059, 
                        RBP = 37134107, 
                        reg = 'inver1')

chr36_inver1 = invClust(roi = ROI_1_chr36, 
                       wh = 1, 
                       geno = chr36_geno, 
                       annot = chr36_map, 
                       dim = 2)


## inversion plots per SNP cluster on the CHR
plot(chr36_inver1)

# 
# ## chr36 sv size 
# 
# chr36_map %>% 
#   slice(301:400) %>% 
#   summarize(start = first(position),
#             end = last(position)) %>% 
#   mutate(sv_size = end-start, 
#          sv_size_mb = sv_size/1000000)
# 
# chr36_inver4$datin$y %>% 
#   as_tibble() %>% 
#   write_csv('chr36_inver4_100SNPS_MDS.csv')
# 
# invGenotypes(chr36_inver4) %>% 
#   as_tibble() %>% 
#   write_csv('chr36_inver4_100SNPS_Inversion_genos.csv')
# 
# 
# # plot chr36  sv --------------------------------------------------------------
# 
# chr36_mds = read_csv('chr36_inver4_100SNPS_MDS.csv')
# chr36_inver_geno = read_csv('chr36_inver4_100SNPS_Inversion_genos.csv')
# chr36_ped = read_table2('AC_New_CHRSET_5.ped', 
#                        col_names = F) %>% 
#   dplyr::select(1:2)
# 
# chr36_inver = bind_cols(chr36_ped, 
#                        chr36_mds, 
#                        chr36_inver_geno)
# 
# 
# shades_of_BigD = c('#0583F2', 
#                             '#F28705',
#                             '#F20530')
#                             
# # View(North_chr1_Big_D_Energy)
# chr36_Big_Plot_Energy = ggplot(data = chr36_inver, 
#                               aes(x = V1, 
#                                   y = V2, 
#                                   col = value))+
#   geom_point()+
#   scale_colour_manual(values = shades_of_BigD)+
#   # facet_grid(.~label)+
#   labs(x = 'MDS axis 1', 
#        y = 'MDS axis 2', 
#        col = 'Inversion genotypes', 
#        title = 'chr36 11.2 Mb')+
#   theme(panel.grid = element_blank(), 
#         axis.title = element_text(size = 14), 
#         axis.text.y = element_text(size = 12), 
#         axis.text.x = element_text(size = 10, 
#                                    angle = 45, 
#                                    hjust = 1, 
#                                    vjust = 1),
#         legend.title = element_text(size = 14), 
#         legend.text = element_text(size = 12), 
#         strip.background = element_rect(fill = 'white'),
#         strip.text = element_text(size = 12, 
#                                   face = 'bold'))
# 
# chr36_Big_Plot_Energy
# 
# ggsave('chr36_putate_sv.tiff', 
#        plot = chr36_Big_Plot_Energy, 
#        dpi = 'retina', 
#        units = 'cm', 
#        width = 20, 
#        height = 15)
















