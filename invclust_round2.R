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


