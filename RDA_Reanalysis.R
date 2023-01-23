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

## sdm predictors dataset
sdm_env = read_csv('sdmpredictors_All_Populations_04.01.2022.csv') %>% 
  dplyr::rename(Population = Name) %>% 
  na.omit()

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


## cleaning the sdm predictors data set
sdm_env_data = left_join(sdm_env, 
                         LatLong, 
                         by = 'Population')
sdm_env_data = left_join(genotype_data, 
                         sdm_env_data, 
                         by = c('Population',
                                'Lat',
                                'Long')) %>%
  # dplyr::select(1:5) %>% 
# View()
  # na.omit() %>% 
  dplyr::select(Population,
                IID,
                Lat,
                Long,
                starts_with('env')) %>% 
  dplyr::rename(icecover = 5, 
                temp_mean = 6, 
                chloro_mean = 7, 
                dissox_mean = 8, 
                iron_mean = 9, 
                phosphate_mean = 10, 
                nitrate_mean = 11, 
                primprod_mean = 12, 
                salinity_mean = 13, 
                silicate_mean = 14) %>% 
  na.omit()




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

## sdmpredictors snp matrix

sdm_SNPs = inner_join(sdm_env_data, 
           genotype_data) %>% 
  dplyr::select(matches('AX.')) %>% 
  as_tibble()


# 
# sdm_SNPs = left_join(genotype_data, 
#           sdm_env_data, 
#           by = c('Population',
#                  'Lat',
#                  'Long')) %>% 
#   dplyr::select(matches('AX.')) %>% 
#   as_tibble()


(sum(is.na(sdm_SNPs))/13123)*100


## impute the NA's to the most common genotype at the locus
sdm_SNPs = apply(sdm_SNPs ,
                     2,
                     function(x) replace(x, 
                                         is.na(x),
                                         as.numeric(names(which.max(table(x))))))
sdm_SNPs = sdm_SNPs %>%
  as_tibble()


# sdm partial rda ---------------------------------------------------------
sdm_partial_RDA = rda(sdm_SNPs ~ icecover + dissox_mean + primprod_mean + Condition(sdm_env_data$Lat), 
                          data = sdm_env_data, 
                          scale = T)

RsquareAdj(sdm_partial_RDA)
summary(eigenvals(sdm_partial_RDA,
                  model = "constrained"))

vif.cca(sdm_partial_RDA)
# sdm_sig = anova.cca(sdm_partial_RDA) 

sdm_partial_sum = summary(sdm_partial_RDA)


# sdm partial rda outliers ------------------------------------------------

## write data to files.
sdm_partial_sum$species %>%
  as_tibble() %>%
  write_csv('AC_sdm_Partial_RDA_All_Pops_SNPS_23.01.2023.csv')

## RDA scores for each individual
sdm_partial_sum$sites %>%
  as_tibble() %>%
  write_csv('AC_sdm_Partial_RDA_All_Pops_INDIVIDUALS_23.01.2023.csv')

## RDA scores for the variables used as predictors
sdm_partial_sum$biplot %>%
  as_tibble() %>%
  write_csv('AC_sdm_Partial_RDA_All_Pops_BIPLOTVARS_23.01.2023.csv')

## All three axes were significant!
## Need to pull snps from all three axes
sdm_partial_RDA_scores = scores(sdm_partial_RDA,
                                    choices = c(1:3),
                                    display = 'species')

## Pulling out the outlier loci for each significant axis
## Outliers are 3 standard deviations from the axis mean

## The first axis explains almost 90% variation
## only using the outliers on the first axis 
## BUT NEED SECOND AXIS FOR THE BIPLOT
sdm_partial_RDA_outliers_axis1 = outliers(sdm_partial_RDA_scores[,1],3)
sdm_partial_RDA_outliers_axis2 = outliers(sdm_partial_RDA_scores[,2],3)
# sdm_partial_RDA_outliers_axis3 = outliers(sdm_partial_RDA_scores[,3],3)

## Creating a cleaned data frame for outiers on each axis
sdm_partial_RDA_out_axis1 = cbind.data.frame(rep(1,
                                                     times = length(sdm_partial_RDA_outliers_axis1)),
                                                 names(sdm_partial_RDA_outliers_axis1),
                                                 unname(sdm_partial_RDA_outliers_axis1))%>%
  as_tibble() %>%
  dplyr::rename(Axis = 1,
                SNP = 2,
                RDA_score = 3)
sdm_partial_RDA_out_axis2 = cbind.data.frame(rep(2,
                                                     times = length(sdm_partial_RDA_outliers_axis2)),
                                                 names(sdm_partial_RDA_outliers_axis2),
                                                 unname(sdm_partial_RDA_outliers_axis2)) %>%
  as_tibble() %>%
  dplyr::rename(Axis = 1,
                SNP = 2,
                RDA_score = 3)
# 
# sdm_partial_RDA_out_axis3 = cbind.data.frame(rep(3,
#                                                      times = length(sdm_partial_RDA_outliers_axis3)),
#                                                  names(sdm_partial_RDA_outliers_axis3),
#                                                  unname(sdm_partial_RDA_outliers_axis3)) %>%
#   as_tibble() %>%
#   dplyr::rename(Axis = 1,
#                 SNP = 2,
#                 RDA_score = 3) 

## Full cleaned data frame with all outliers
## used all three axes in the RDA
sdm_partial_RDA_out_total = bind_rows(sdm_partial_RDA_out_axis1,
                                      sdm_partial_RDA_out_axis2)

# sdm_partial_RDA_out_total = bind_rows(sdm_partial_RDA_out_axis1,
#                                           sdm_partial_RDA_out_axis2,
#                                           sdm_partial_RDA_out_axis3)


## the rda scores for all snps, not just the outliers
sdm_all_snps = sdm_partial_RDA_scores[,1] %>% 
  as.data.frame()

## data frame for the normy snps
sdm_partial_RDA_normy = cbind.data.frame(rep(0,
                                                 times = length(sdm_all_snps)),
                                             row.names(sdm_all_snps),
                                             unname(sdm_all_snps))

sdm_partial_RDA_normy = sdm_partial_RDA_normy %>%
  as_tibble() %>%
  dplyr::rename(Axis = 1,
                SNP = 2,
                RDA_score_axis1 = 3)

# View(sdm_partial_RDA_normy)
## If this doesn't work, you might need to load the data.table R package
## this pulls out the outlier snps from the full snp data frame
## we only want the normal nonoutlier snps
sdm_partial_RDA_normy = sdm_partial_RDA_normy[!sdm_partial_RDA_normy$SNP %in% sdm_partial_RDA_out_axis1$SNP,]
sdm_partial_RDA_normy = sdm_partial_RDA_normy[!sdm_partial_RDA_normy$SNP %in% sdm_partial_RDA_out_axis2$SNP,]
# sdm_partial_RDA_normy = sdm_partial_RDA_normy[!sdm_partial_RDA_normy$SNP %in% sdm_partial_RDA_out_axis3$SNP,]

write_csv(sdm_partial_RDA_normy,
          'AC_sdm_partial_RDA_Associations_Normy_SNPs_23.01.2023.csv')


## get the predictor variables associated with each outlier locus
# sdm_partial_RDA_out_total = as.data.frame(sdm_partial_RDA_out_total)
sdm_SNPs = as.data.frame(sdm_SNPs)

sdm_env = sdm_env_data %>% 
  dplyr::select(Lat, 
                Long, 
                icecover, 
                dissox_mean, 
                primprod_mean) %>% 
  as.data.frame()

sdm_partial_RDA_outliers = sdm_partial_RDA_out_total %>% 
  as.data.frame()

nam = sdm_partial_RDA_outliers[1:492, 2]
# nam = RDA_out[1:109,2]
out_snps = sdm_SNPs[,nam]
outlier_correlations = apply(sdm_env, 
                             2, 
                             function(x)cor(x, 
                                            out_snps))
out_snp_cor = as_tibble(outlier_correlations)

sdm_partial_RDA_outliers= bind_cols(sdm_partial_RDA_outliers,
                                        out_snp_cor) %>%
  as_tibble()

# View(sdm_partial_RDA_out_total)

## check for duplicated outliers across axes
# length(sdm_partial_RDA_out_total$SNP[duplicated(sdm_partial_RDA_out_total$SNP)])

## for loop to get the predctor and correlation coefficient 
## for each variable in the RDA
for (i in 1:length(sdm_partial_RDA_outliers$SNP)) {
  bar <- sdm_partial_RDA_outliers[i,]
  sdm_partial_RDA_outliers[i,9] <- names(which.max(abs(bar[6:8]))) # gives the variable
  sdm_partial_RDA_outliers[i,10] <- max(abs(bar[6:8]))              # gives the correlation
}

## clean up the new dataframe
sdm_cand_snps = sdm_partial_RDA_outliers %>% 
  rename(predictor = ...9, 
         correlation = ...10)

# View(sdm_cand_snps)

## save the data
write_csv(sdm_cand_snps, 
          'AC_sdm_partial_RDA_Association_Outlier_SNPs_23.01.2023.csv')


##


# sdm fix snp labels ------------------------------------------------------


## SNP formats don't match and need to be updated to 
## correspond to the map file
map = read_tsv('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.map', 
               col_names = c('Chromosome', 
                             'SNP', 
                             'Genetic_pos', 
                             'Physical_pos'))
sdm_partial_outs = read_csv('AC_sdm_partial_RDA_Association_Outlier_SNPs_23.01.2023.csv')

sdm_normy_snps = read_csv('AC_sdm_partial_RDA_Associations_Normy_SNPs_23.01.2023.csv')

## This gets rid the format the snps are in after the RDA
## We need to line up the SNP names to the map file
## The first set of two operations arefor the outlier snps
sdm_partial_outs$SNP = gsub("AX.",
                                "AX-",
                                sdm_partial_outs$SNP)

sdm_partial_outs$SNP = gsub("_.*",
                                "",
                                sdm_partial_outs$SNP)

## The next set of operations are for the non-outlier snps
sdm_normy_snps$SNP = gsub("AX.",
                              "AX-",
                              sdm_normy_snps$SNP)

sdm_normy_snps$SNP = gsub("_.*",
                              "",
                              sdm_normy_snps$SNP)

## This gets the map file data for the outlier snps
sdm_outs_map = map[map$SNP %in% sdm_partial_outs$SNP,]
## This merges the map file with the outlier snps
sdm_partial_outs = merge(sdm_outs_map,
                             sdm_partial_outs,
                             by.x = 'SNP',
                             by.y = 'SNP') %>%
  as_tibble()

## Now we're going to do the same thing with the non-outlier snps
normy_map = map[map$SNP %in% sdm_normy_snps$SNP,]
sdm_normy_snps = merge(normy_map,
                           sdm_normy_snps,
                           by.x = 'SNP',
                           by.y = 'SNP') %>%
  as_tibble()

## Write out the rda scores for all of the map data
write_csv(sdm_normy_snps,
          'AC_sdm_partial_RDA_Normysnp_data_23.01.2023.csv')

## This gets us the rda scores for all of the snps used

# map_all_snp = map %>%
#   dplyr::select(SNP)
# sdm_Full_scores = as.data.frame(cbind(SNP = rownames(sdm_partial_RDA_scores),
#                                           sdm_partial_RDA_scores)) %>%
#   as_tibble()
# sdm_Full_scores = bind_cols(map_all_snp,
#                                 sdm_Full_scores)
# 
# sdm_out_snps_rdascores = merge(sdm_partial_outs,
#                                    sdm_Full_scores,
#                                    by.x = 'SNP',
#                                    by.y = 'SNP') %>%
#   as_tibble()
# 
# write_csv(sdm_out_snps_rdascores,
#           'AC_sdm_partial_RDA_outlier_data_08.06.2022.csv')
#


##Chromosome 9 (LG6.2) and 16 (LG13) have >20 outliers which is different 
## than what I found before
sdm_partial_outs %>% 
  arrange(Chromosome) %>% 
  group_by(Chromosome) %>% 
  summarise(n = n()) %>% 
  # filter(n > 20) %>% 
  View()


# sdm partial biplot ------------------------------------------------------



# sdm normal rda ----------------------------------------------------------

sdm_RDA = rda(sdm_SNPs ~ icecover + dissox_mean + primprod_mean, 
                      data = sdm_env_data, 
                      scale = T)

RsquareAdj(sdm_RDA)
summary(eigenvals(sdm_RDA,
                  model = "constrained"))

vif.cca(sdm_RDA)
# sdm_sig = anova.cca(sdm_partial_RDA) 

sdm_sum = summary(sdm_RDA)


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


##Chromosome 9 (LG6.2) and 16 (LG13) have >20 outliers which is different 
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
bioclim_sum = summary(bioclim_RDA)

# bioclim  rda outliers --------------------------------------------

## write data to files.
bioclim_sum$species %>%
  as_tibble() %>%
  write_csv('AC_Bioclim_RDA_All_Pops_SNPS_09.01.2023.csv')

## RDA scores for each individual
bioclim_sum$sites %>%
  as_tibble() %>%
  write_csv('AC_Bioclim_RDA_All_Pops_INDIVIDUALS_09.01.2023.csv')

## RDA scores for the variables used as predictors
bioclim_sum$biplot %>%
  as_tibble() %>%
  write_csv('AC_Bioclim_RDA_All_Pops_BIPLOTVARS_09.01.2023.csv')

## All three axes were significant!
## Need to pull snps from all three axes
bioclim_RDA_scores = scores(bioclim_RDA,
                                    choices = c(1:3),
                                    display = 'species')

## Pulling out the outlier loci for each significant axis
## Outliers are 3 standard deviations from the axis mean

## The first axis explains almost 90% variation
## only using the outliers on the first axis
bioclim_RDA_outliers = outliers(bioclim_RDA_scores[,1],3)
# bioclim_partial_RDA_outliers_axis2 = outliers(bioclim_partial_RDA_scores[,2],3)
# bioclim_partial_RDA_outliers_axis3 = outliers(bioclim_partial_RDA_scores[,3],3)

## Creating a cleaned data frame for outiers on each axis
bioclim_RDA_out = cbind.data.frame(rep(1,
                                                     times = length(bioclim_RDA_outliers)),
                                                 names(bioclim_RDA_outliers),
                                                 unname(bioclim_RDA_outliers))%>%
  as_tibble() %>%
  dplyr::rename(Axis = 1,
                SNP = 2,
                RDA_score = 3)


## the rda scores for all snps, not just the outliers
bioclim_all_snps = bioclim_RDA_scores[,1] %>% 
  as.data.frame()

## data frame for the normy snps
bioclim_RDA_normy = cbind.data.frame(rep(0,
                                                 times = length(bioclim_all_snps)),
                                             row.names(bioclim_all_snps),
                                             unname(bioclim_all_snps))

bioclim_RDA_normy = bioclim_RDA_normy %>%
  as_tibble() %>%
  dplyr::rename(Axis = 1,
                SNP = 2,
                RDA_score_axis1 = 3)


## If this doesn't work, you might need to load the data.table R package
## this pulls out the outlier snps from the full snp data frame
## we only want the normal nonoutlier snps
bioclim_RDA_normy = bioclim_RDA_normy[!bioclim_RDA_normy$SNP %in% bioclim_RDA_out$SNP,]
# bioclim_partial_RDA_normy = bioclim_partial_RDA_normy[!bioclim_partial_RDA_normy$SNP %in% bioclim_partial_RDA_out_axis2$SNP,]
# bioclim_partial_RDA_normy = bioclim_partial_RDA_normy[!bioclim_partial_RDA_normy$SNP %in% bioclim_partial_RDA_out_axis3$SNP,]

write_csv(bioclim_RDA_normy,
          'AC_RDA_Associations_Normy_SNPs_09.01.2023.csv')


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

bioclim_RDA_outliers = bioclim_RDA_out %>% 
  as.data.frame()

nam = bioclim_RDA_outliers[1:139, 2]
# nam = RDA_out[1:109,2]
out_snps = bioclim_SNPS[,nam]
outlier_correlations = apply(bioclim_env, 
                             2, 
                             function(x)cor(x, 
                                            out_snps))
out_snp_cor = as_tibble(outlier_correlations)

bioclim_RDA_outliers= bind_cols(bioclim_RDA_outliers,
                                        out_snp_cor) %>%
  as_tibble()

# View(bioclim_partial_RDA_out_total)

## check for duplicated outliers across axes
# length(bioclim_partial_RDA_out_total$SNP[duplicated(bioclim_partial_RDA_out_total$SNP)])

## for loop to get the predctor and correlation coefficient 
## for each variable in the RDA

for (i in 1:length(bioclim_RDA_outliers$SNP)) {
  bar <- bioclim_RDA_outliers[i,]
  bioclim_RDA_outliers[i,9] <- names(which.max(abs(bar[6:8]))) # gives the variable
  bioclim_RDA_outliers[i,10] <- max(abs(bar[6:8]))              # gives the correlation
}

## clean up the new dataframe
bioclim_cand_snps = bioclim_RDA_outliers %>% 
  rename(predictor = ...9, 
         correlation = ...10)

# View(bioclim_cand_snps)

## save the data
write_csv(bioclim_cand_snps, 
          'AC_RDA_Association_Outlier_SNPs_09.01.2023.csv')

# Fix SNP labels ----------------------------------------------------------

## SNP formats don't match and need to be updated to 
## correspond to the map file
map = read_tsv('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.map', 
               col_names = c('Chromosome', 
                             'SNP', 
                             'Genetic_pos', 
                             'Physical_pos'))
bioclim_outs = read_csv('AC_RDA_Association_Outlier_SNPs_09.01.2023.csv')

bioclim_normy_snps = read_csv('AC_RDA_Associations_Normy_SNPs_09.01.2023.csv')

## This gets rid the format the snps are in after the RDA
## We need to line up the SNP names to the map file
## The first set of two operations arefor the outlier snps
bioclim_outs$SNP = gsub("AX.",
                                "AX-",
                                bioclim_outs$SNP)

bioclim_outs$SNP = gsub("_.*",
                                "",
                                bioclim_outs$SNP)

## The next set of operations are for the non-outlier snps
bioclim_normy_snps$SNP = gsub("AX.",
                              "AX-",
                              bioclim_normy_snps$SNP)

bioclim_normy_snps$SNP = gsub("_.*",
                              "",
                              bioclim_normy_snps$SNP)

## This gets the map file data for the outlier snps
bioclim_outs_map = map[map$SNP %in% bioclim_outs$SNP,]
## This merges the map file with the outlier snps
bioclim_outs = merge(bioclim_outs_map,
                             bioclim_outs,
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
          'AC_bioclim_RDA_Normysnp_data_09.01.2023.csv')

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
bioclim_outs %>% 
  arrange(Chromosome) %>% 
  group_by(Chromosome) %>% 
  summarise(n = n()) %>% 
  # filter(n > 20) %>% 
  View()




# Bioclim partial RDA biplot ----------------------------------------------
theme_set(theme_bw())

bioclim_rda_snps = read_csv('AC_Bioclim_Partial_RDA_All_Pops_SNPS_09.01.2023.csv')
bioclim_rda_indivs = read_csv('AC_Bioclim_Partial_RDA_All_Pops_INDIVIDUALS_09.01.2023.csv')
bioclim_rda_vars = read_csv('AC_Bioclim_Partial_RDA_All_Pops_BIPLOTVARS_09.01.2023.csv')
# bio_clim_data = read_csv('bioclim_data_mito_nuclear_charr.csv')
bioclim_rda_env_data = read_csv('bioclim_data_AC_AllPops.csv') %>% 
  dplyr::select(Long, 
                Lat, 
                bio1, 
                bio3, 
                bio4)
bioclim_rda_outliers = read_csv('AC_partial_RDA_Association_Outlier_SNPs_09.01.2023.csv')
bioclim_rda_normy_snps = read_csv('AC_partial_RDA_Associations_Normy_SNPs_09.01.2023.csv')

# rda_outliers %>% 
#   filter(Axis != 3) %>% 
#   View()
bioclim_col.pred = rownames(bioclim_partial_RDA$CCA$v) %>% 
  as_tibble() %>% 
  dplyr::rename(SNP = value)

bioclim_col.pred$SNP = gsub("AX.",
                            "AX-",
                            bioclim_col.pred$SNP)

bioclim_col.pred$SNP = gsub("_.*",
                            "",
                            bioclim_col.pred$SNP)

bioclim_sel = bioclim_rda_outliers$SNP
bioclim_env = bioclim_rda_outliers$predictor
bioclim_env_col = bioclim_rda_outliers$predictor
bioclim_env_col[bioclim_env_col == 'bio1'] = '#663F8C'
bioclim_env_col[bioclim_env_col == 'bio3'] = '#D9965B'
bioclim_env_col[bioclim_env_col == 'bio4'] = '#BF4B54'

bioclim_sel_tib = bioclim_sel %>% 
  as_tibble()
bioclim_env_tib = bioclim_env %>% 
  as_tibble()
bioclim_env_col_tib = bioclim_env_col %>% 
  as_tibble()

bioclim_SNP_cols = bind_cols(bioclim_sel_tib, 
                             bioclim_env_tib, 
                             bioclim_env_col_tib) %>% 
  dplyr::rename(SNP = 1, 
                var = 2, 
                col = 3)

bioclim_rda_out_test = bind_cols(bioclim_rda_outliers, 
                                 bioclim_SNP_cols)


bioclim_RDA_biplot = ggplot() +
  ## this geom point is for the zoomed out version
  geom_point(data = bioclim_rda_snps,
             aes(x = RDA1,
                 y = RDA2),
             col = '#A6A6A6', 
             size = 2)+
  geom_point(data = bioclim_rda_out_test,
             aes(x = RDA_score_axis1,
                 y = RDA_score_axis2,
                 col = var),
             size = 2)+
  # geom_point(data = rda_outliers, 
  #            aes(x = RDA_score_axis1, 
  #                y = RDA_score_axis2), 
  #            col = '#0AB33A', 
  #            size = 2)+
  ## These two things are only for the zoomed out rda to look
  ## at population level variation related to isotopic variation
# geom_point(data = bioclim_rda_env_data,
#            aes(x = RDA1,
#                y = RDA2,
#                col = Latitude),
#            size = 2)+
# scale_color_viridis(option = 'magma')+
# scale_color_manual(values = c('#663F8C',
#                               '#D9965B',
#                               '#BF4B54'))+
scale_color_manual(values = c('#A6036D',
                              '#F2D6A2',
                              '#8DF2CD'))+
  ## THis is for the zoomed out RDA to show var in pops
  geom_segment(aes(xend = bioclim_RDA$CCA$biplot[,1],
                   yend = bioclim_RDA$CCA$biplot[,2],
                   x=0,
                   y=0),
               colour="black",
               size=1,
               linetype=1,
               arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x = 1.0*bioclim_RDA$CCA$biplot[,1],
                y = 1.2*bioclim_RDA$CCA$biplot[,2],
                label = colnames(bioclim_rda_env_data[,3:5]))) +
  
  labs(x = 'RDA 1 (83.7% variance explained)',
       y = 'RDA 2 (13.1% variance explained)', 
       title = 'A)')+
  theme(#legend.position = "none", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 14), 
    axis.ticks = element_line(size = 1), 
    plot.title = element_text(size = 15, 
                              hjust = 0), 
    # legend.title = element_text(size = 13),
    # legend.text = element_text(size = 12),
    legend.position = 'none'
  )

bioclim_RDA_biplot


# Bioclim partial manhattan plot ------------------------------------------


