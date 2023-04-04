##################################################################
## PCAdmix file processing
##
## Matt Brachmann (PhDMattyB)
##
## 2021-03-25
##
##################################################################

## Load the data manipulation work horse
library(tidyverse)
library(superheat)
library(data.table)
library(reshape2)
library(viridis)


# setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/')

setwd('~/Bradbury_Postdoc/AC_SV_MS_Data/RFmix/')

# Plink data --------------------------------------------------------------


map = read_tsv('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.map', 
               col_names = F) %>% 
  rename(`#CHROMOSOME` = X1, 
         MARKERID = X2, 
         GENETIC_DIST = X3, 
         PHYSICAL_DIST = X4)

data = read_table('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.ped', 
                  col_names = F)

data %>% 
  select(IndividualID)
##
# chr convert -------------------------------------------------------------

# Chr_convert = function(data){
#   data = mutate(.data = data,
#                 AC_CHR = as.factor(case_when(
#                   CHROMOSOME == 'NC_036838.1' ~ '1',
#                   CHROMOSOME == 'NC_036839.1' ~ '2',
#                   CHROMOSOME == 'NC_036840.1' ~ '3',
#                   CHROMOSOME == 'NC_036841.1' ~ '4',
#                   CHROMOSOME == 'NC_036842.1' ~ '5',
#                   CHROMOSOME == 'NC_036843.1' ~ '6',
#                   CHROMOSOME == 'NC_036844.1' ~ '7',
#                   CHROMOSOME == 'NC_036845.1' ~ '8',
#                   CHROMOSOME == 'NC_036846.1' ~ '9',
#                   CHROMOSOME == 'NC_036847.1' ~ '10',
#                   CHROMOSOME == 'NC_036848.1' ~ '11',
#                   CHROMOSOME == 'NC_036849.1' ~ '12',
#                   CHROMOSOME == 'NC_036850.1' ~ '13',
#                   CHROMOSOME == 'NC_036851.1' ~ '14',
#                   CHROMOSOME == 'NC_036852.1' ~ '15',
#                   CHROMOSOME == 'NC_036853.1' ~ '16',
#                   CHROMOSOME == 'NC_036854.1' ~ '17',
#                   CHROMOSOME == 'NC_036855.1' ~ '18',
#                   CHROMOSOME == 'NC_036856.1' ~ '19',
#                   CHROMOSOME == 'NC_036857.1' ~ '20',
#                   CHROMOSOME == 'NC_036858.1' ~ '21',
#                   CHROMOSOME == 'NC_036859.1' ~ '22',
#                   CHROMOSOME == 'NC_036860.1' ~ '23',
#                   CHROMOSOME == 'NC_036861.1' ~ '24',
#                   CHROMOSOME == 'NC_036862.1' ~ '25',
#                   CHROMOSOME == 'NC_036863.1' ~ '26',
#                   CHROMOSOME == 'NC_036864.1' ~ '27',
#                   CHROMOSOME == 'NC_036865.1' ~ '28',
#                   CHROMOSOME == 'NC_036866.1' ~ '29',
#                   CHROMOSOME == 'NC_036867.1' ~ '30',
#                   CHROMOSOME == 'NC_036868.1' ~ '31',
#                   CHROMOSOME == 'NC_036869.1' ~ '32',
#                   CHROMOSOME == 'NC_036870.1' ~ '33',
#                   CHROMOSOME == 'NC_036871.1' ~ '34',
#                   CHROMOSOME == 'NC_036872.1' ~ '35',
#                   CHROMOSOME == 'NC_036873.1' ~ '36',
#                   CHROMOSOME == 'NC_036874.1' ~ '37',
#                   CHROMOSOME == 'NC_036875.1' ~ '38',
#                   CHROMOSOME == 'NC_036876.1' ~ '39',
#                   CHROMOSOME > 'NC_036876.1' ~ '40')))
#   return(data)  
# }
# 
# map = Chr_convert(map) %>% 
#   select(AC_CHR, 
#          MARKERID, 
#          GENETIC_DIST, 
#          PHYSICAL_DIST) %>% 
#   rename(`#Chromosome` = AC_CHR)




# Make --keep files for admixed populations -------------------------------

ped_test = read_table2('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.ped', 
                       col_names = F)

ped_ids = read_table2('Charr_Poly_All_Fixed_coords_maf05_geno95.fam', 
                      col_names = F) %>%
  dplyr::select(X1,
                X2)

set.seed(666)

ATL = ped_test %>% 
  filter(X1 %in% c('VTG', 
                            'FLJ', 
                            'GAL', 
                            'MJO',
                            'SVI', 
                            'THI', 
                            'VAT', 
                            'LBN', 
                            'DUG'))

Iceland_sample = ATL %>% 
  filter(X1 %in% c('FLJ',
                            'GAL', 
                            'MJO', 
                            'SVI', 
                            'THI', 
                            'VAT')) %>% 
  sample_n(size = 30, 
           replace = F)

ATL = ATL %>% 
  filter(X1 %in% c('VTG', 
                            'LBN', 
                            'DUG')) %>% 
  bind_rows(Iceland_sample) %>% 
  sample_n(size = 20, 
           replace = F)

ACD = data %>% 
  filter(X1 %in% c('GNL', 
                            'PEN')) %>% 
  sample_n(size = 20, 
           replace = F)

ARC = data %>% 
  filter(X1 %in% c('BTR', 
                            'HAB', 
                            'IPI', 
                            'SGL', 
                            'AUP', 
                            'BRT',
                            'FMI')) %>% 
  sample_n(size = 20, 
           replace = F)



dim(ATL)
dim(ACD)
dim(ARC)

## Filter the data to be the admixed populations
admixed = anti_join(data, 
                    ATL, 
                    by = 'X1')

admixed = anti_join(admixed, 
                    ACD, 
                    by = 'X1')
admixed = anti_join(admixed, 
                    ARC, 
                    by = 'X1')

admixed %>% dplyr::select(1) %>% 
  distinct() %>% 
  View()
dim(admixed)

ATL %>% 
  dplyr::select(X1, 
                X2) %>% 
  write_tsv('ATL_ref_keep.txt', 
            col_names = F)

ARC %>% 
  dplyr::select(X1, 
                X2) %>% 
  write_tsv('ARC_ref_keep.txt', 
            col_names = F)

ACD %>% 
  dplyr::select(X1, 
                X2) %>% 
  write_tsv('ACD_ref_keep.txt', 
            col_names = F)
admixed %>% 
 
  dplyr::select(X1, 
                X2) %>% 
  filter(! X1 %in% c('SVI', 
                     'THI', 
                     'GAL', 
                     'FLJ')) 

ind_sample = data %>% 
  dplyr::select(X2)
pop_sample = data %>% 
  dplyr::select(X1)

sample_map = bind_cols(ind_sample, 
                       pop_sample) 

ACD %>% 
  dplyr::select(X1) %>% 
  distinct() %>% 
  arrange()
# %>% 
#   write_tsv('rfmix_sample_map.txt')
data = mutate(.data = sample_map,
              pop = as.factor(case_when(
                                X1 == 'ANA' ~ 'admixed',
                                X1 == 'IKI' ~ 'admixed2',
                                X1 == 'IKA' ~ 'admixed',
                                X1 == 'KIN' ~ 'admixed',
                                X1 == 'PAL' ~ 'admixed',
                                X1 == 'STV' ~ 'admixed',
                                X1 == 'R103' ~ 'admixed',
                                X1 == 'KIY' ~ 'admixed',
                                X1 == 'R78' ~ 'admixed',
                                X1 == 'GDL' ~ 'admixed',
                                X1 == 'BLD' ~ 'admixed',
                                X1 == 'KOM' ~ 'admixed',
                                X1 == 'FRD' ~ 'admixed',
                                X1 == 'SWA' ~ 'admixed',
                                X1 == 'PUT' ~ 'admixed',
                                X1 == 'AVA' ~ 'admixed',
                                X1 == 'R104' ~ 'admixed',
                                X1 == 'KOG' ~ 'admixed',
                                X1 == 'MBB' ~ 'admixed',
                                X1 == 'R110' ~ 'admixed',
                                X1 == 'FRS' ~ 'admixed',
                                X1 == 'PAN' ~ 'admixed',
                                X1 == 'NAK' ~ 'admixed',
                                X1 == 'FRN' ~ 'admixed',
                                X1 == 'R109' ~ 'admixed',
                                X1 == 'NAC' ~ 'admixed',
                                X1 == 'KAM' ~ 'admixed',
                                X1 == 'ENG' ~ 'admixed',
                                X1 == 'TOR' ~ 'admixed',
                                X1 == 'KAN' ~ 'admixed',
                                X1 == 'IKL' ~ 'admixed',
                                X1 == 'MCC' ~ 'admixed',
                                X1 == 'UNH' ~ 'admixed',
                                X1 == 'NOR' ~ 'admixed',
                                X1 == 'BRG' ~ 'admixed',
                                X1 == 'REI' ~ 'admixed',
                                X1 == 'R105' ~ 'admixed',
                                X1 == 'PBP' ~ 'admixed',
                                X1 == 'IGL' ~ 'admixed',
                                X1 == 'SGL' ~ 'ARC',
                                X1 == 'BRT' ~ 'ARC',
                                X1 == 'IPI' ~ 'ARC',
                                X1 == 'AUP' ~ 'ARC',
                                X1 == 'HAB' ~ 'ARC',
                                X1 == 'FMI' ~ 'ARC',
                                X1 == 'BTR' ~ 'ARC',
                                X1 == 'LBN' ~ 'ATL',
                                X1 == 'VTG' ~ 'ATL',
                                X1 == 'DUG' ~ 'ATL',
                                X1 == 'VAT' ~ 'ATL', 
                                X1 == 'MJO' ~ 'ATL',
                                X1 == 'GNL' ~ 'ACD', 
                                X1 == 'PEN' ~ 'ACD')))

sample_map = data %>% 
  dplyr::select(X2, 
                pop)

write.table(x = sample_map, 
            file = 'rfmix_sample_map.txt', 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')


chr = map %>% 
  dplyr::select(1)
phys = map %>% 
  dplyr::select(4)
genet = map %>% 
  dplyr::select(3)

genetic_map = bind_cols(chr, phys, genet) 

write.table(x = genetic_map, 
            file = 'rfmix_genetic_map.txt', 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

##
# Separate by reference and admixed populations ---------------------------

## Making the ATL, ACD, ARC reference populations
## for the three glacial lineages
ATL = data %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'FLJ', 
                            'GAL', 
                            'MJO',
                            'SVI', 
                            'THI', 
                            'VAT', 
                            'LBN', 
                            'DUG'))

Iceland_sample = ATL %>% 
  filter(`#FamilyID` %in% c('FLJ',
                            'GAL', 
                            'MJO', 
                            'SVI', 
                            'THI', 
                            'VAT')) %>% 
  sample_n(size = 30, 
           replace = F)

ATL = ATL %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'LBN', 
                            'DUG')) %>% 
  bind_rows(Iceland_sample)


ACD = data %>% 
  filter(`#FamilyID` %in% c('GNL', 
                            'PEN'))

ARC = data %>% 
  filter(`#FamilyID` %in% c('BTR', 
                            'HAB', 
                            'IPI', 
                            'SGL', 
                            'AUP', 
                            'BRT',
                            'FMI')) %>% 
  sample_n(size = 100, 
           replace = F)



dim(ATL)
dim(ACD)
dim(ARC)

## Filter the data to be the admixed populations
admixed = anti_join(data, 
          ATL, 
          by = '#FamilyID')

admixed = anti_join(admixed, 
                    ACD, 
                    by = '#FamilyID')
admixed = anti_join(admixed, 
                    ARC, 
                    by = '#FamilyID')

## CHeck against our list to double check
admixed %>% 
  dplyr::select(`#FamilyID`) %>% 
  distinct() %>% 
  arrange(`#FamilyID`) %>% 
  View()


# Karas population subset ------------------------------------------------
ATL = data %>% 
  filter(`#FamilyID` == 'LBN')
ACD = data %>% 
  filter(`#FamilyID` == 'PEN')
ARC = data %>% 
  filter(`#FamilyID` == 'IPI')

dim(ATL)
dim(ACD)
dim(ARC)


# Ancestral subset round 2 ------------------------------------------------
## having an issue with there being to many reference
## individuals than there are markers on some chromosomes
## randomly subsample down to 50 individuals per reference
## population

ATL = data %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'FLJ', 
                            'GAL', 
                            'MJO',
                            'SVI', 
                            'THI', 
                            'VAT', 
                            'LBN', 
                            'DUG'))

Iceland_sample = ATL %>% 
  filter(`#FamilyID` %in% c('FLJ',
                            'GAL', 
                            'MJO', 
                            'SVI', 
                            'THI', 
                            'VAT')) %>% 
  sample_n(size = 30, 
           replace = F)

ATL = ATL %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'LBN', 
                            'DUG')) %>% 
  bind_rows(Iceland_sample) %>% 
  sample_n(size = 50, 
           replace = F)


ACD = data %>% 
  filter(`#FamilyID` %in% c('GNL', 
                            'PEN'))

ARC = data %>% 
  filter(`#FamilyID` %in% c('BTR', 
                            'HAB', 
                            'IPI', 
                            'SGL', 
                            'AUP', 
                            'BRT',
                            'FMI')) %>% 
  sample_n(size = 50, 
           replace = F)



dim(ATL)
dim(ACD)
dim(ARC)

## Filter the data to be the admixed populations
admixed = anti_join(data, 
                    ATL, 
                    by = '#FamilyID')

admixed = anti_join(admixed, 
                    ACD, 
                    by = '#FamilyID')
admixed = anti_join(admixed, 
                    ARC, 
                    by = '#FamilyID')

## CHeck against our list to double check
admixed %>% 
  select(`#FamilyID`) %>% 
  distinct() %>% 
  arrange(`#FamilyID`) %>% 
  View()

# Ancestral subset round 3 - 40 individuals ------------------------------------------------
## having an issue with there being to many reference
## individuals than there are markers on some chromosomes
## randomly subsample down to 50 individuals per reference
## population

ATL = data %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'FLJ', 
                            'GAL', 
                            'MJO',
                            'SVI', 
                            'THI', 
                            'VAT', 
                            'LBN', 
                            'DUG'))

Iceland_sample = ATL %>% 
  filter(`#FamilyID` %in% c('FLJ',
                            'GAL', 
                            'MJO', 
                            'SVI', 
                            'THI', 
                            'VAT')) %>% 
  sample_n(size = 30, 
           replace = F)

ATL = ATL %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'LBN', 
                            'DUG')) %>% 
  bind_rows(Iceland_sample) %>% 
  sample_n(size = 40, 
           replace = F)

ACD = data %>% 
  filter(`#FamilyID` %in% c('GNL', 
                            'PEN')) %>% 
  sample_n(size = 40, 
           replace = F)

ARC = data %>% 
  filter(`#FamilyID` %in% c('BTR', 
                            'HAB', 
                            'IPI', 
                            'SGL', 
                            'AUP', 
                            'BRT',
                            'FMI')) %>% 
  sample_n(size = 40, 
           replace = F)



dim(ATL)
dim(ACD)
dim(ARC)

## Filter the data to be the admixed populations
# admixed = anti_join(data, 
#                     ATL, 
#                     by = '#FamilyID')
# 
# admixed = anti_join(admixed, 
#                     ACD, 
#                     by = '#FamilyID')
# admixed = anti_join(admixed, 
#                     ARC, 
#                     by = '#FamilyID')

## CHeck against our list to double check
# admixed %>% 
#   select(`#FamilyID`) %>% 
#   distinct() %>% 
#   arrange(`#FamilyID`) %>% 
#   View()


# Ancestral subset round 4 - 30 individuals ------------------------------------------------
## having an issue with there being to many reference
## individuals than there are markers on some chromosomes
## randomly subsample down to 50 individuals per reference
## population

ATL = data %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'FLJ', 
                            'GAL', 
                            'MJO',
                            'SVI', 
                            'THI', 
                            'VAT', 
                            'LBN', 
                            'DUG'))

Iceland_sample = ATL %>% 
  filter(`#FamilyID` %in% c('FLJ',
                            'GAL', 
                            'MJO', 
                            'SVI', 
                            'THI', 
                            'VAT')) %>% 
  sample_n(size = 30, 
           replace = F)

ATL = ATL %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'LBN', 
                            'DUG')) %>% 
  bind_rows(Iceland_sample) %>% 
  sample_n(size = 30, 
           replace = F)

ACD = data %>% 
  filter(`#FamilyID` %in% c('GNL', 
                            'PEN')) %>% 
  sample_n(size = 30, 
           replace = F)

ARC = data %>% 
  filter(`#FamilyID` %in% c('BTR', 
                            'HAB', 
                            'IPI', 
                            'SGL', 
                            'AUP', 
                            'BRT',
                            'FMI')) %>% 
  sample_n(size = 30, 
           replace = F)



dim(ATL)
dim(ACD)
dim(ARC)



# Ancestral subset round 5 - 20 individuals ------------------------------------------------
## having an issue with there being to many reference
## individuals than there are markers on some chromosomes
## randomly subsample down to 50 individuals per reference
## population
set.seed(666)
ATL = data %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'FLJ', 
                            'GAL', 
                            'MJO',
                            'SVI', 
                            'THI', 
                            'VAT', 
                            'LBN', 
                            'DUG'))

Iceland_sample = ATL %>% 
  filter(`#FamilyID` %in% c('FLJ',
                            'GAL', 
                            'MJO', 
                            'SVI', 
                            'THI', 
                            'VAT')) %>% 
  sample_n(size = 30, 
           replace = F)

ATL = ATL %>% 
  filter(`#FamilyID` %in% c('VTG', 
                            'LBN', 
                            'DUG')) %>% 
  bind_rows(Iceland_sample) %>% 
  sample_n(size = 20, 
           replace = F)

ACD = data %>% 
  filter(`#FamilyID` %in% c('GNL', 
                            'PEN')) %>% 
  sample_n(size = 20, 
           replace = F)

ARC = data %>% 
  filter(`#FamilyID` %in% c('BTR', 
                            'HAB', 
                            'IPI', 
                            'SGL', 
                            'AUP', 
                            'BRT',
                            'FMI')) %>% 
  sample_n(size = 20, 
           replace = F)



dim(ATL)
dim(ACD)
dim(ARC)

## Filter the data to be the admixed populations
admixed = anti_join(data, 
                    ATL, 
                    by = '#FamilyID')

admixed = anti_join(admixed, 
                    ACD, 
                    by = '#FamilyID')
admixed = anti_join(admixed, 
                    ARC, 
                    by = '#FamilyID')
dim(admixed)


##

# Ancestral no greenland 20 refs ------------------------------------------
# Ancestral subset round 5 - 20 individuals ------------------------------------------------
## having an issue with there being to many reference
## individuals than there are markers on some chromosomes
## randomly subsample down to 50 individuals per reference
## population
set.seed(666)
ATL = data %>% 
  filter(`#FamilyID` %in% c(
    # 'VTG', 
                            'FLJ', 
                            'GAL', 
                            'MJO',
                            'SVI', 
                            'THI', 
                            'VAT', 
                            'LBN', 
                            'DUG'))

Iceland_sample = ATL %>% 
  filter(`#FamilyID` %in% c('FLJ',
                            'GAL', 
                            'MJO', 
                            'SVI', 
                            'THI', 
                            'VAT')) %>% 
  sample_n(size = 30, 
           replace = F)

ATL = ATL %>% 
  filter(`#FamilyID` %in% c(
    # 'VTG', 
                            'LBN', 
                            'DUG')) %>% 
  bind_rows(Iceland_sample) %>% 
  sample_n(size = 20, 
           replace = F)

ACD = data %>% 
  filter(`#FamilyID` %in% c('GNL', 
                            'PEN')) %>% 
  sample_n(size = 20, 
           replace = F)

ARC = data %>% 
  filter(`#FamilyID` %in% c('BTR', 
                            'HAB', 
                            'IPI', 
                            'SGL', 
                            'AUP', 
                            'BRT',
                            'FMI')) %>% 
  sample_n(size = 20, 
           replace = F)



dim(ATL)
dim(ACD)
dim(ARC)

## Filter the data to be the admixed populations
admixed = anti_join(data, 
                    ATL, 
                    by = '#FamilyID')

admixed = anti_join(admixed, 
                    ACD, 
                    by = '#FamilyID')
admixed = anti_join(admixed, 
                    ARC, 
                    by = '#FamilyID')
dim(admixed)


##


# Write tsv files to put back into plink ----------------------------------

## write these plink files to the pcadmix folder
ATL %>% 
  write_tsv('ATL_ref_20.ped')
ACD %>% 
  write_tsv('ACD_ref_20.ped')
ARC %>% 
  write_tsv('ARC_ref_20.ped')
admixed %>%
  write_tsv('Admixed_populations.ped', 
            col_names = F)

admixed %>% 
  dplyr::select(1:15) %>% 
  slice(1:6) %>% 
  View()

map %>%
  write_tsv('ATL_ref_20.map')
map %>%
  write_tsv('ACD_ref_20.map')
map %>%
  write_tsv('ARC_ref_20.map')
map %>%
  write_tsv('Admixed_populations.map')

# Create Individual ID column to match vit --------------------------------
admixed = read_tsv('Admixed_populations.ped')

adm = admixed %>% 
  dplyr::select(IndividualID) %>% 
  separate(col = IndividualID, 
           into = c('pop', 'id'), 
           sep = '_')

pop = adm %>% 
  dplyr::select(pop) %>% 
  as.data.frame()
id = adm %>% 
  dplyr::select(id) %>% 
  as.data.frame()

adm1 = paste0(pop[1:nrow(pop),], 
              '_A_', 
              id[1:nrow(id),]) %>% 
  as_tibble() %>% 
  rename(IndividualID = value)

adm2 = paste0(pop[1:nrow(pop),], 
              '_B_', 
              id[1:nrow(id),]) %>% 
  as_tibble() %>% 
  rename(IndividualID = value)

adm_id = bind_rows(adm1, 
                   adm2) %>%
  separate(col = IndividualID, 
           into = c('pop', 
                    'id1', 
                    'id2'),
           sep = '_') %>% 
  arrange(pop, id2) %>% 
  unite(IndividualID, 
        c('pop', 
          'id1', 
          'id2'), 
        sep = '_')
##
# PCAdmix results ---------------------------------------------------------
setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_20individuals_ref/')

read_tsv('Admixed_populations_chr-1.map', 
         col_names = c('Chr', 
                       'Marker', 
                       'Genetic_dist', 
                       'Physical_dist'))

read_table2('PCAdmix_20individuals_results_chr-1.vit.txt')

read_table2('PCAdmix_20individuals_results_chr-1.fbk.txt')
fread('PCAdmix_20individuals_results_chr-1.fbk.txt')

read_table2('PCAdmix_20individuals_results_chr-1.ia.txt')
fread('PCAdmix_20individuals_results_chr-1.ia.txt')

read_table2('PCAdmix_20individuals_results_chr-1.markers.txt')

read_table2('PCAdmix_20individuals_results_chr-1.pca.txt')


# Results Kara run individuals --------------------------------------------------

setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_kararun_ref/')

vit_file = fread('PCAdmix_kararun_chr-11.vit.txt')%>% 
  as_tibble()

dim(vit_file)

vit_matrix = as.matrix(vit_file[,2:length(vit_file)])
superheat(vit_matrix)
## 0 = ATL
## 1 = ARC
## 2 = ACD


##



# Results 50 individuals --------------------------------------------------

vit_file = fread('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAadmix_50individuals_ref/PCAdmix_50individuals_results_chr11.vit.txt') %>% 
  as_tibble()
# Results 20NoGreenland ---------------------------------------------------
vit_file = fread('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_ref_20NoGreenland/PCAdmix_20NoGreenland_results_chr-11.vit.txt')%>% 
  as_tibble()


# Results 20 individuals --------------------------------------------------
setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_20individuals_ref/')

## contribution of each ancestor population to each window
fread('PCAdmix_20individuals_results_chr-11.fbk.txt')

## I honestly don't know what this file is.....
fread('PCAdmix_20individuals_results_chr-11.ia.txt')

## markers used in the windows
# fread('PCAdmix_20individuals_results_chr-11.markers.txt') %>% 
#   as_tibble()

Markers = read_table2('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_20individuals_ref/PCAdmix_20individuals_results_chr-1.markers.txt', 
            col_names = F)

## pca loadings for the chromosome
pca = fread('PCAdmix_20individuals_results_chr-11.pca.txt') %>% 
  as_tibble() 

# ggplot(data = pca, 
#        aes(x = V2, 
#            y = V3, 
#            col = V1))+
#   geom_point()

vit_file = fread('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_20individuals_ref/PCAdmix_20individuals_results_Chr-11.vit.txt')%>% 
  as_tibble()

# dim(vit_file)

# vit_matrix = as.matrix(vit_file[,2:length(vit_file)])
# superheat(vit_matrix)
## 0 = ATL
## 1 = ARC
## 2 = ACD

vit_data = bind_cols(vit_file, 
          adm_id) %>% 
  select(IndividualID, 
         everything())
##
# average the data per individual--------------------------------------------------------

avg_vit_data = stats::filter(vit_file[,2:length(vit_file)], 
         rep(1/2, 
             2), 
             sides = 1) %>%
    as.data.frame() %>% 
    as_tibble() %>% 
    filter(row_number() %% 2 == 0)

identifiers = admixed %>%
  dplyr::select(`#FamilyID`,
         IndividualID)

clean_vit_data = bind_cols(identifiers,
          avg_vit_data)

write_csv(clean_vit_data, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_ref_20NoGreenland/Chr-11_Averaged_PCAdmix_data.csv')

# averaged Heatmap plots -----------------------------------------------------------
# clean_vit_data = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_20individuals_ref/Chr-39_Averaged_PCAdmix_data.csv')
clean_vit_data = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_ref_20NoGreenland/Chr-11_Averaged_PCAdmix_data.csv')

clean_melted_vit = melt(clean_vit_data, 
     id.vars = c('#FamilyID', 'IndividualID')) %>% 
  as_tibble()  

# melted_vit$IndividualID = as.factor(melted_vit$IndividualID)
# melted_vit$IndividualID = as.numeric(as.character(melted_vit$IndividualID))
theme_set(theme_bw())

glacial_cols = c('#4E94BF', 
                 '#0D5A8F',
                 '#4E458C',
                 '#B050E0',
                 '#F23545')

glacial_cols = c('#F23545',
                 '#B050E0',
                 '#4E458C',
                 '#0D5A8F',
                 '#4E94BF')

clean_melted_vit$value = as.character(clean_melted_vit$value)

pcadmix = clean_melted_vit %>% 
  ggplot(aes(x = variable,
         y = IndividualID,
         col = value, 
         fill = value)) +
  geom_bar(stat = 'identity', 
           width = 1)+
  scale_fill_manual(values = glacial_cols)+
  scale_color_manual(values = glacial_cols)+
  # scale_fill_viridis(viridis(n=5))+
  # scale_color_viridis(viridis(n=5))+
  labs(x = 'Window', 
       y = 'Individuals', 
       colour = 'Glacial lineage')+
  theme(axis.text.y = element_blank(),
        axis.title = element_text(size = 14),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 12, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1), 
        legend.title = element_blank(),
        # legend.position = 'none'
        )

pcadmix_popsplit = clean_melted_vit %>% 
  ggplot(aes(x = variable, 
             y = IndividualID, 
             col = value, 
             fill = value)) +
  geom_bar(stat = 'identity', 
           width = 1)+
  scale_fill_viridis(viridis(n=5))+
  scale_color_viridis(viridis(n=5))+
  facet_grid(~ `#FamilyID`)+
  labs(x = 'Window', 
       y = 'Individuals')+
  theme(axis.text = element_blank(), 
        axis.title = element_text(size = 14),
        panel.grid = element_blank(), 
        strip.text = element_text(face = 'bold'),
        strip.background = element_rect(fill = 'white'), 
        legend.title = element_blank())

pcadmix
pcadmix_popsplit
# dir.create('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/MKB_per_chr_plots_NoGreenland')
# 

## DOUBLE CHECK THE FOLDER YOU'RE SAVING THE RESULTS IN!!!!!!

ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/MKB_per_chr_plots/Chr-11_pcadmix_20refs_glacialcols.tiff',
      plot = pcadmix,
      dpi = 'retina',
      units = 'cm', 
      width = 18, 
      height = 11)

ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/MKB_per_chr_plots/Chr-39_pcadmix_20refs_popsplit.tiff',
       plot = pcadmix_popsplit,
       dpi = 'retina',
       units = 'cm',
       width = 45,
       height = 10)


# ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/MKB_per_chr_plots_NoGreenland/Chr-11_pcadmix_20refs.tiff', 
#        plot = pcadmix, 
#        dpi = 'retina', 
#        units = 'cm')
# 
# ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/MKB_per_chr_plots_NoGreenland/Chr-11_pcadmix_20refs_popsplit.tiff', 
#        plot = pcadmix_popsplit, 
#        dpi = 'retina', 
#        units = 'cm', 
#        width = 45, 
#        height = 10)

# ggsave('PCAdmix_KaraRun_Chr-11.tiff', 
#        plot = pcadmix, 
#        dpi = 'retina', 
#        units = 'cm')
# 
# ggsave('PCAdmix_PopSplit_KaraRun_Chr-11.tiff', 
#        plot = pcadmix_popsplit, 
#        dpi = 'retina', 
#        units = 'cm', 
#        width = 45, 
#        height = 10)

#
# heatmap no averaging (OG results) ---------------------------------------
theme_set(theme_bw())
library(reshape2)
library(viridis)
melted_vit = melt(vit_data, 
                  id.vars = c('V1', 'IndividualID')) %>% 
  as_tibble()  

# melted_vit$IndividualID = as.factor(melted_vit$IndividualID)
# melted_vit$IndividualID = as.numeric(as.character(melted_vit$IndividualID))

OG_pcadmix_chr11 = melted_vit %>% 
  ggplot(aes(x = variable, 
             y = IndividualID, 
             col = value, 
             fill = value)) +
  geom_bar(stat = 'identity', 
           width = 1)+
  scale_fill_viridis(viridis(n=5))+
  scale_color_viridis(viridis(n=5))+
  labs(x = 'Window', 
       y = 'Individuals', 
       colour = 'Glacial lineage')+
  theme(axis.text.y = element_blank(),
        axis.title = element_text(size = 14),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 12, 
                                   angle = 45), 
        legend.title = element_blank())

OG_pcadmix_chr11_pops = melted_vit %>% 
  ggplot(aes(x = variable, 
             y = IndividualID, 
             col = value, 
             fill = value)) +
  geom_bar(stat = 'identity', 
           width = 1)+
  scale_fill_viridis(viridis(n=5))+
  scale_color_viridis(viridis(n=5))+
  facet_grid(~ V1)+
  labs(x = 'Window', 
       y = 'Individuals')+
  theme(axis.text = element_blank(), 
        axis.title = element_text(size = 14),
        panel.grid = element_blank(), 
        strip.text = element_text(face = 'bold'),
        strip.background = element_rect(fill = 'white'), 
        legend.title = element_blank())


ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Figures/PCAdmix_chr11_OG_results.tiff',
       plot = OG_pcadmix_chr11,
       dpi = 'retina',
       units = 'cm', 
       width = 25, 
       height = 15)

ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Figures/PCAdmix_chr11_OG_results_per_pop.tiff',
       plot = OG_pcadmix_chr11_pops,
       dpi = 'retina',
       units = 'cm',
       width = 80,
       height = 10)
##
# match lostruct outliers -------------------------------------------------
# chr 11

clean_vit_data %>% 
  select(`#FamilyID`, 
         IndividualID, 
         V6, 
         V7, 
         V13, 
         V14) %>% 
  write_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/PCAdmix_outlier_ancestry.csv')

##
# Resuts 30 individuals ---------------------------------------------------
setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_30individuals_ref/')
vit_file = fread('PCAdmix_30individuals_results_chr-11.vit.txt')%>% 
  as_tibble()

dim(vit_file)

vit_matrix = as.matrix(vit_file[,2:length(vit_file)])
superheat(vit_matrix)
## 0 = ATL
## 1 = ARC
## 2 = ACD

# Resuts 40 individuals ---------------------------------------------------
setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_40individuals_ref/')
vit_file = fread('PCAdmix_40individuals_results_chr-11.vit.txt')%>% 
  as_tibble()

dim(vit_file)

vit_matrix = as.matrix(vit_file[,2:length(vit_file)])
superheat(vit_matrix)
## 0 = ATL
## 1 = ARC
## 2 = ACD

# Resuts 50 individuals ---------------------------------------------------
setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAadmix_50individuals_ref/')
vit_file = fread('PCAdmix_50individuals_results_chr11.vit.txt')%>% 
  as_tibble()

dim(vit_file)

vit_matrix = as.matrix(vit_file[,2:length(vit_file)])
superheat(vit_matrix)
## 0 = ATL
## 1 = ARC
## 2 = ACD
# PCAdmix plots -----------------------------------------------------------

setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/PCAdmix/PCAdmix_w5_adm_Lab/')


ac08_data = read_tsv('NC_036848_PCAdmix.pca.txt')


ac08_vit = read_csv('NC_036848_PCAdmix.vit.csv')
ac08_markers = read_tsv('NC_036848_PCAdmix.markers.txt')


read.table("NC_036848_PCAdmix.vit.csv", 
           row.names = 'V1')
