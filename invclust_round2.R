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
plot(chr23_inver3)


## chr23 sv size 

chr23_map %>% 
  slice(301:400) %>% 
  summarize(start = first(position),
            end = last(position)) %>% 
  mutate(sv_size = end-start, 
         sv_size_mb = sv_size/1000000)

chr23_inver4$datin$y %>% 
  as_tibble() %>% 
  write_csv('chr23_inver4_100SNPS_MDS.csv')

invGenotypes(chr23_inver4) %>% 
  as_tibble() %>% 
  write_csv('chr23_inver4_100SNPS_Inversion_genos.csv')


# plot chr23  sv --------------------------------------------------------------

chr23_mds = read_csv('chr23_inver4_100SNPS_MDS.csv')
chr23_inver_geno = read_csv('chr23_inver4_100SNPS_Inversion_genos.csv')
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
       title = 'chr23 11.2 Mb')+
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



