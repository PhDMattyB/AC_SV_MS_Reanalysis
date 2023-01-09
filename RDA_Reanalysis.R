##############################
## AC RDA analysis
##
## Matt Brachmann (PhDMattyB)
##
## 2023-01-09
##
##############################

setwd('~/AC_SV_MS_Data/RDA_Analysis/')

library(tidyverse)
library(raster)
library(MASS)
library(vegan)
library(psych)
library(sdmpredictors)
# library(leaflet)
# library(ggnewscale)
library(patchwork)

theme_set(theme_bw())

## Vignette for the RDA association analysis
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html

## Stolen from Forester et al., 2018 tutorial on RDA analysis
## This function defines the outlier loci
outliers = function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]]
  # locus names in these tails
}



# Load Data ---------------------------------------------------------------


## lat long
LatLong = read_csv('SampleSites_Coords_1June2020.csv') 

## bioclim data
bioclim_env_data = read_csv('bioclim_data_AC_AllPops.csv')

## genotype data
genotype_data = read.delim('Charr_Poly_All_Fixed_coords_maf05_geno95_RecodeA.raw', 
                           sep = "", 
                           stringsAsFactors = F) %>% 
  as_tibble() %>% 
  dplyr::rename(Population = FID)



# Clean Data --------------------------------------------------------------

## Need to add the lat long data to the genotype data
genotype_data = left_join(genotype_data, 
                          LatLong, 
                          by = 'Population') %>% 
  dplyr::select(Population,
                IID, 
                Lat, 
                Long, 
                starts_with('AX.'))
## Need the lat long for all the bioclimatic variables
bioclim_env_data = left_join(genotype_data,
                             bioclim_env_data,
                             by = c('Lat', 
                                    'Long')) %>%
  # na.omit() %>% 
  dplyr::select(Population,
                IID,
                Lat,
                Long,
                starts_with('bio'))


# bioclim_env_data %>% 
#   distinct(Population, .keep_all = T) %>% 
#   View()

## we don't have bioclimatic variables for 3/57 populations

bioclim_data_naomit = bioclim_env_data %>% 
  na.omit()

bioclim_full_data = inner_join(bioclim_data_naomit, 
                               genotype_data)


# Create SNP matrix -------------------------------------------------------

## Need to create a data frame with just the SNPs we want to look at
bioclim_SNPS = bioclim_full_data %>%
  dplyr::select(matches("AX.")) %>%
  as_tibble(bioclim_full_data)

dim(bioclim_SNPS)
## The realignment lost us about 4000 SNPS, great

## Check for the number of NA's
(sum(is.na(bioclim_SNPS))/13123)*100


## impute the NA's to the most common genotype at the locus
bioclim_SNPS = apply(bioclim_SNPS ,
                     2,
                     function(x) replace(x, 
                                         is.na(x),
                                         as.numeric(names(which.max(table(x))))))
bioclim_SNPS = bioclim_SNPS %>%
  as_tibble()


# Bioclim partial RDA -----------------------------------------------------

bioclim_partial_RDA = rda(bioclim_SNPS ~ bio1 + bio3 + bio4 + Condition(bioclim_data_naomit$Lat), 
                          data = bioclim_data_naomit, 
                          scale = T)

RsquareAdj(bioclim_partial_RDA)
summary(eigenvals(bioclim_partial_RDA,
                  model = "constrained"))

vif.cca(bioclim_partial_RDA)
# bioclim_sig = anova.cca(bioclim_partial_RDA) 

bioclim_partial_sum = summary(bioclim_partial_RDA)


# bioclim partial rda outliers --------------------------------------------

## write data to files.
bioclim_partial_sum$species %>%
  as_tibble() %>%
  write_csv('AC_Bioclim_Partial_RDA_All_Pops_SNPS_09.01.2023.csv')

## RDA scores for each individual
bioclim_partial_sum$sites %>%
  as_tibble() %>%
  write_csv('AC_Bioclim_Partial_RDA_All_Pops_INDIVIDUALS_09.01.2023.csv')

## RDA scores for the variables used as predictors
bioclim_partial_sum$biplot %>%
  as_tibble() %>%
  write_csv('AC_Bioclim_Partial_RDA_All_Pops_BIPLOTVARS_09.01.2023.csv')

## All three axes were significant!
## Need to pull snps from all three axes
bioclim_partial_RDA_scores = scores(bioclim_partial_RDA,
                                    choices = c(1:3),
                                    display = 'species')

## Pulling out the outlier loci for each significant axis
## Outliers are 3 standard deviations from the axis mean

## The first axis explains almost 90% variation
## only using the outliers on the first axis
bioclim_partial_RDA_outliers_axis1 = outliers(bioclim_partial_RDA_scores[,1],3)
# bioclim_partial_RDA_outliers_axis2 = outliers(bioclim_partial_RDA_scores[,2],3)
# bioclim_partial_RDA_outliers_axis3 = outliers(bioclim_partial_RDA_scores[,3],3)

## Creating a cleaned data frame for outiers on each axis
bioclim_partial_RDA_out_axis1 = cbind.data.frame(rep(1,
                                                     times = length(bioclim_partial_RDA_outliers_axis1)),
                                                 names(bioclim_partial_RDA_outliers_axis1),
                                                 unname(bioclim_partial_RDA_outliers_axis1))%>%
  as_tibble() %>%
  dplyr::rename(Axis = 1,
                SNP = 2,
                RDA_score = 3)
# bioclim_partial_RDA_out_axis2 = cbind.data.frame(rep(2,
#                                                      times = length(bioclim_partial_RDA_outliers_axis2)),
#                                                  names(bioclim_partial_RDA_outliers_axis2),
#                                                  unname(bioclim_partial_RDA_outliers_axis2)) %>%
#   as_tibble() %>%
#   dplyr::rename(Axis = 1,
#                 SNP = 2,
#                 RDA_score = 3) 
# 
# bioclim_partial_RDA_out_axis3 = cbind.data.frame(rep(3,
#                                                      times = length(bioclim_partial_RDA_outliers_axis3)),
#                                                  names(bioclim_partial_RDA_outliers_axis3),
#                                                  unname(bioclim_partial_RDA_outliers_axis3)) %>%
#   as_tibble() %>%
#   dplyr::rename(Axis = 1,
#                 SNP = 2,
#                 RDA_score = 3) 

## Full cleaned data frame with all outliers
## used all three axes in the RDA

# bioclim_partial_RDA_out_total = bind_rows(bioclim_partial_RDA_out_axis1, 
#                                           bioclim_partial_RDA_out_axis2, 
#                                           bioclim_partial_RDA_out_axis3)


## the rda scores for all snps, not just the outliers
bioclim_all_snps = bioclim_partial_RDA_scores[,1] %>% 
  as.data.frame()

## data frame for the normy snps
bioclim_partial_RDA_normy = cbind.data.frame(rep(0,
                                                 times = length(bioclim_all_snps)),
                                             row.names(bioclim_all_snps),
                                             unname(bioclim_all_snps))

bioclim_partial_RDA_normy = bioclim_partial_RDA_normy %>%
  as_tibble() %>%
  dplyr::rename(Axis = 1,
                SNP = 2,
                RDA_score_axis1 = 3)


## If this doesn't work, you might need to load the data.table R package
## this pulls out the outlier snps from the full snp data frame
## we only want the normal nonoutlier snps
bioclim_partial_RDA_normy = bioclim_partial_RDA_normy[!bioclim_partial_RDA_normy$SNP %in% bioclim_partial_RDA_out_axis1$SNP,]
# bioclim_partial_RDA_normy = bioclim_partial_RDA_normy[!bioclim_partial_RDA_normy$SNP %in% bioclim_partial_RDA_out_axis2$SNP,]
# bioclim_partial_RDA_normy = bioclim_partial_RDA_normy[!bioclim_partial_RDA_normy$SNP %in% bioclim_partial_RDA_out_axis3$SNP,]

write_csv(bioclim_partial_RDA_normy,
          'AC_partial_RDA_Associations_Normy_SNPs_09.01.2023.csv')


## get the predictor variables associated with each outlier locus
# bioclim_partial_RDA_out_total = as.data.frame(bioclim_partial_RDA_out_total)
bioclim_SNPS = as.data.frame(bioclim_SNPS)

bioclim_env = bioclim_full_data %>% 
  dplyr::select(Lat, 
                Long, 
                bio1, 
                bio3, 
                bio4) %>% 
  as.data.frame()

bioclim_partial_RDA_outliers = bioclim_partial_RDA_out_axis1 %>% 
  as.data.frame()

nam = bioclim_partial_RDA_outliers[1:140, 2]
# nam = RDA_out[1:109,2]
out_snps = bioclim_SNPS[,nam]
outlier_correlations = apply(bioclim_env, 
                             2, 
                             function(x)cor(x, 
                                            out_snps))
out_snp_cor = as_tibble(outlier_correlations)

bioclim_partial_RDA_outliers= bind_cols(bioclim_partial_RDA_outliers,
                                          out_snp_cor) %>%
  as_tibble()

# View(bioclim_partial_RDA_out_total)

## check for duplicated outliers across axes
# length(bioclim_partial_RDA_out_total$SNP[duplicated(bioclim_partial_RDA_out_total$SNP)])

## for loop to get the predctor and correlation coefficient 
## for each variable in the RDA
for (i in 1:length(bioclim_partial_RDA_outliers$SNP)) {
  bar <- bioclim_partial_RDA_outliers[i,]
  bioclim_partial_RDA_outliers[i,9] <- names(which.max(abs(bar[6:8]))) # gives the variable
  bioclim_partial_RDA_outliers[i,10] <- max(abs(bar[6:8]))              # gives the correlation
}

## clean up the new dataframe
bioclim_cand_snps = bioclim_partial_RDA_outliers %>% 
  rename(predictor = ...9, 
         correlation = ...10)

# View(bioclim_cand_snps)

## save the data
write_csv(bioclim_cand_snps, 
          'AC_partial_RDA_Association_Outlier_SNPs_09.01.2023.csv')

# Fix SNP labels ----------------------------------------------------------

## SNP formats don't match and need to be updated to 
## correspond to the map file
map = read_tsv('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.map', 
               col_names = c('Chromosome', 
                             'SNP', 
                             'Genetic_pos', 
                             'Physical_pos'))
bioclim_partial_outs = read_csv('AC_partial_RDA_Association_Outlier_SNPs_09.01.2023.csv')

bioclim_normy_snps = read_csv('AC_partial_RDA_Associations_Normy_SNPs_09.01.2023.csv')

## This gets rid the format the snps are in after the RDA
## We need to line up the SNP names to the map file
## The first set of two operations arefor the outlier snps
bioclim_partial_outs$SNP = gsub("AX.",
                            "AX-",
                            bioclim_partial_outs$SNP)

bioclim_partial_outs$SNP = gsub("_.*",
                            "",
                            bioclim_partial_outs$SNP)

## The next set of operations are for the non-outlier snps
bioclim_normy_snps$SNP = gsub("AX.",
                              "AX-",
                              bioclim_normy_snps$SNP)

bioclim_normy_snps$SNP = gsub("_.*",
                              "",
                              bioclim_normy_snps$SNP)

## This gets the map file data for the outlier snps
bioclim_outs_map = map[map$SNP %in% bioclim_partial_outs$SNP,]
## This merges the map file with the outlier snps
bioclim_partial_outs = merge(bioclim_outs_map,
                         bioclim_partial_outs,
                         by.x = 'SNP',
                         by.y = 'SNP') %>%
  as_tibble()

## Now we're going to do the same thing with the non-outlier snps
normy_map = map[map$SNP %in% bioclim_normy_snps$SNP,]
bioclim_normy_snps = merge(normy_map,
                           bioclim_normy_snps,
                           by.x = 'SNP',
                           by.y = 'SNP') %>%
  as_tibble()

## Write out the rda scores for all of the map data
write_csv(bioclim_normy_snps,
          'AC_bioclim_partial_RDA_Normysnp_data_09.01.2023.csv')

## This gets us the rda scores for all of the snps used

# map_all_snp = map %>%
#   dplyr::select(SNP)
# bioclim_Full_scores = as.data.frame(cbind(SNP = rownames(bioclim_partial_RDA_scores),
#                                           bioclim_partial_RDA_scores)) %>%
#   as_tibble()
# bioclim_Full_scores = bind_cols(map_all_snp,
#                                 bioclim_Full_scores)
# 
# bioclim_out_snps_rdascores = merge(bioclim_partial_outs,
#                                    bioclim_Full_scores,
#                                    by.x = 'SNP',
#                                    by.y = 'SNP') %>%
#   as_tibble()
# 
# write_csv(bioclim_out_snps_rdascores,
#           'AC_bioclim_partial_RDA_outlier_data_08.06.2022.csv')
#


##Chromosome 9 and 16 have >20 outliers which is different 
## than what I found before
bioclim_partial_outs %>% 
  arrange(Chromosome) %>% 
  group_by(Chromosome) %>% 
  summarise(n = n()) %>% 
  # filter(n > 20) %>% 
  View()


# bioclim normal rda ------------------------------------------------------


## VIF gets better without the partial RDA
## overall r-squared is higher without the condition
## The model has a much larger effect. 
## generally both RDA's show similar patterns
## I can plot that to demonstrate that
bioclim_RDA = rda(bioclim_SNPS ~ bio1 + bio3 + bio4, 
                          data = bioclim_data_naomit, 
                          scale = T)

RsquareAdj(bioclim_RDA)
summary(eigenvals(bioclim_RDA,
                  model = "constrained"))

vif.cca(bioclim_RDA)
# bioclim_sig = anova.cca(bioclim_RDA) 
bioclim_rda_sum = summary(bioclim_RDA)


