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
                       qvalues)


# Graph admixture ---------------------------------------------------------


melted_dwata = melt(Clean_data, 
                   id.vars = c('X1', 
                               'X2')) %>% 
  as_tibble()

## bw is the king of plots
theme_set(theme_bw())

## need a colour palette
test_col = c( '#4E9EBF',
              '#4E458C',
              '#F23545', 
              '#F29F05')


admixture = ggplot(data = melted_dwata, 
                      aes(x = X2,
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
        #                            size = 6,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

admixture
