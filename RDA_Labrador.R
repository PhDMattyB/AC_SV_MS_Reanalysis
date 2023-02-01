## RDA analysis for the Newfoundland and Labrador populations only

setwd('~/Bradbury_Postdoc/AC_SV_MS_Data/RDA_Analysis/')

library(tidyverse)
library(raster)
library(MASS)
library(vegan)
library(psych)
library(sdmpredictors)
# library(leaflet)
# library(ggnewscale)
library(patchwork)
library(ggnewscale)

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



# Load data ---------------------------------------------------------------


## lat long
LatLong = read_csv('SampleSites_Coords_1June2020.csv') 

## bioclim data
bioclim_env_data = read_csv('bioclim_data_AC_AllPops.csv')

## genotype data
genotype_data = read.delim('Charr_Poly_All_Fixed_coords_maf05_geno95_envmatch_Lab.raw', 
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
(sum(is.na(bioclim_SNPS))/12596)*100


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


(sum(is.na(sdm_SNPs))/12596)*100


## impute the NA's to the most common genotype at the locus
sdm_SNPs = apply(sdm_SNPs ,
                 2,
                 function(x) replace(x, 
                                     is.na(x),
                                     as.numeric(names(which.max(table(x))))))
sdm_SNPs = sdm_SNPs %>%
  as_tibble()


# sdm partial rda ---------------------------------------------------------
# sdm_env_data %>% 
#   group_by(Population) %>% 
#   summarize(number = n()) %>% 
#   View()

sdm_partial_RDA = rda(sdm_SNPs ~ icecover + dissox_mean + primprod_mean + Condition(sdm_env_data$Lat), 
                      data = sdm_env_data, 
                      scale = T)

RsquareAdj(sdm_partial_RDA)
summary(eigenvals(sdm_partial_RDA,
                  model = "constrained"))

vif.cca(sdm_partial_RDA)
# sdm_sig = anova.cca(sdm_partial_RDA) 

sdm_partial_sum = summary(sdm_partial_RDA)


# sdm lab outliers --------------------------------------------------------


## write data to files.
sdm_partial_sum$species %>%
  as_tibble() %>%
  write_csv('LAB_AC_sdm_Partial_RDA_All_Pops_SNPS_01.02.2023.csv')

## RDA scores for each individual
sdm_partial_sum$sites %>%
  as_tibble() %>%
  write_csv('LAB_AC_sdm_Partial_RDA_All_Pops_INDIVIDUALS_01.02.2023.csv')

## RDA scores for the variables used as predictors
sdm_partial_sum$biplot %>%
  as_tibble() %>%
  write_csv('LAB_AC_sdm_Partial_RDA_All_Pops_BIPLOTVARS_01.02.2023.csv')

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

## Full cleaned data frame with all outliers
## used all three axes in the RDA
sdm_partial_RDA_out_total = bind_rows(sdm_partial_RDA_out_axis1,
                                      sdm_partial_RDA_out_axis2)

## the rda scores for all snps, not just the outliers
## data frame for the normy snps

sdm_partial_RDA_normy = sdm_partial_RDA_scores[,1:2] %>% 
  as.data.frame() %>% 
  as_tibble(rownames = 'SNP') %>% 
  mutate(Axis = 0) %>% 
  dplyr::select(Axis, 
                SNP, 
                RDA1, 
                RDA2) %>% 
  rename(RDA_score_axis1 = RDA1, 
         RDA_score_axis2 = RDA2)

# View(sdm_partial_RDA_normy)
## If this doesn't work, you might need to load the data.table R package
## this pulls out the outlier snps from the full snp data frame
## we only want the normal nonoutlier snps
sdm_partial_RDA_normy = sdm_partial_RDA_normy[!sdm_partial_RDA_normy$SNP %in% sdm_partial_RDA_out_axis1$SNP,]
sdm_partial_RDA_normy = sdm_partial_RDA_normy[!sdm_partial_RDA_normy$SNP %in% sdm_partial_RDA_out_axis2$SNP,]
# sdm_partial_RDA_normy = sdm_partial_RDA_normy[!sdm_partial_RDA_normy$SNP %in% sdm_partial_RDA_out_axis3$SNP,]

write_csv(sdm_partial_RDA_normy,
          'LAB_AC_sdm_partial_RDA_Associations_Normy_SNPs_01.02.2023.csv')


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

# sdm_partial_RDA_outliers = sdm_partial_RDA_out_total %>% 
#   as.data.frame()

sdm_partial_RDA_out_axis1 = sdm_partial_RDA_out_total %>%
  filter(Axis == 1) %>% 
  as.data.frame()

nam = sdm_partial_RDA_out_axis1[1:33, 2]
# nam = RDA_out[1:109,2]
out_snps = sdm_SNPs[,nam]
outlier_correlations = apply(sdm_env, 
                             2, 
                             function(x)cor(x, 
                                            out_snps))
out_snp_cor = as_tibble(outlier_correlations)

sdm_partial_RDA_out_axis1= bind_cols(sdm_partial_RDA_out_axis1,
                                     out_snp_cor) %>%
  as_tibble()


sdm_partial_RDA_out_axis2 = sdm_partial_RDA_out_total %>%
  filter(Axis == 2) %>% 
  as.data.frame()

nam = sdm_partial_RDA_out_axis2[1:45, 2]
# nam = RDA_out[1:109,2]
out_snps = sdm_SNPs[,nam]
outlier_correlations = apply(sdm_env, 
                             2, 
                             function(x)cor(x, 
                                            out_snps))
out_snp_cor = as_tibble(outlier_correlations)

sdm_partial_RDA_out_axis2 = bind_cols(sdm_partial_RDA_out_axis2,
                                      out_snp_cor) %>%
  as_tibble()

# View(sdm_partial_RDA_out_total)

## check for duplicated outliers across axes
length(sdm_partial_RDA_out_axis1$SNP[duplicated(sdm_partial_RDA_out_axis2$SNP)])

sdm_partial_RDA_outliers = bind_rows(sdm_partial_RDA_out_axis1, 
                                     sdm_partial_RDA_out_axis2)

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
          'LAB_AC_sdm_partial_RDA_Association_Outlier_SNPs_01.02.2023.csv')

# LAB sdm fix snp labels ------------------------------------------------------


## SNP formats don't match and need to be updated to 
## correspond to the map file
# map = read_tsv('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.map', 
#                col_names = c('Chromosome', 
#                              'SNP', 
#                              'Genetic_pos', 
#                              'Physical_pos'))
map = read_tsv('Lab_map.map', 
               col_names = c('Chromosome', 
                             'SNP', 
                             'Genetic_pos', 
                             'Position'))

sdm_partial_outs = read_csv('LAB_AC_sdm_partial_RDA_Association_Outlier_SNPs_01.02.2023.csv')

sdm_normy_snps = read_csv('LAB_AC_sdm_partial_RDA_Associations_Normy_SNPs_01.02.2023.csv')

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
          'LAB_AC_sdm_partial_RDA_Normysnp_data_01.02.2023.csv')

## This gets us the rda scores for all of the snps used

map_all_snp = map %>%
  dplyr::select(SNP)
sdm_Full_scores = as.data.frame(cbind(SNP = rownames(sdm_partial_RDA_scores),
                                      sdm_partial_RDA_scores)) %>%
  as_tibble()
sdm_Full_scores = bind_cols(map_all_snp,
                            sdm_Full_scores)

sdm_out_snps_rdascores = merge(sdm_partial_outs,
                               sdm_Full_scores,
                               by.x = 'SNP',
                               by.y = 'SNP...1') %>%
  as_tibble()

write_csv(sdm_out_snps_rdascores,
          'LAB_AC_sdm_partial_RDA_outlier_data_01.02.2023.csv')



##Chromosome 9 (LG6.2) and 16 (LG13) have >20 outliers which is different 
## than what I found before
sdm_out_enrich = sdm_partial_outs %>% 
  arrange(Chromosome) %>% 
  group_by(Chromosome) %>% 
  summarise(n = n()) 
# filter(n > 20) %>%
View(sdm_out_enrich)







