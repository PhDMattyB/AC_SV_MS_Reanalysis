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



