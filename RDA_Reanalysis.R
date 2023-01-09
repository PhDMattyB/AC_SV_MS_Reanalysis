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
genotype_data = read.delim('Charr_Poly_All_Fixed_coords_maf05_geno95_envmatch_Lab.raw', 
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
(sum(is.na(bioclim_SNPS))/12596)*100


## impute the NA's to the most common genotype at the locus
bioclim_SNPS = apply(bioclim_SNPS ,
                     2,
                     function(x) replace(x, 
                                         is.na(x),
                                         as.numeric(names(which.max(table(x))))))
bioclim_SNPS = bioclim_SNPS %>%
  as_tibble()


# Bioclim Partial RDA -----------------------------------------------------

bioclim_partial_RDA = rda(bioclim_SNPS ~ bio1 + bio3 + bio4 + Condition(bioclim_env_data$Lat), 
                          data = bioclim_env_data, 
                          scale = T)

RsquareAdj(bioclim_partial_RDA)
summary(eigenvals(bioclim_partial_RDA,
                  model = "constrained"))

vif.cca(bioclim_partial_RDA)
bioclim_sig = anova.cca(bioclim_partial_RDA) 

## VIF gets better without the partial RDA
## overall r-squared is higher without the condition
## The model has a much larger effect. 
## generally both RDA's show similar patterns
## I can plot that to demonstrate that
bioclim_RDA = rda(bioclim_SNPS ~ bio1 + bio3 + bio4, 
                          data = bioclim_env_data, 
                          scale = T)

RsquareAdj(bioclim_RDA)
summary(eigenvals(bioclim_RDA,
                  model = "constrained"))

vif.cca(bioclim_RDA)
bioclim_sig = anova.cca(bioclim_RDA) 


