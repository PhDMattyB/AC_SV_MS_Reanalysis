## Fixing the Arctic charr chromosome naming issue
## the old genome had three names for each of the chromosomes
## the new genome has them ordered 1-39
## we need to use the SNP names to make sure
## things are going to the right spot on the new genome

setwd('~/Bradbury_Postdoc/AC_SV_MS_Data/Charr_SV_MS_data')


new_map= read_tsv('newAC_genome_locations_complete.map', 
               col_names = c('Chromosome', 
                             'SNP', 
                             'Genetic_pos', 
                             'Physical_pos'))
new_map_clean = as.numeric(new_map$Physical_pos) %>% 
  as_tibble() %>% 
  rename(Position = value) %>% 
  bind_cols(., 
            new_map) %>% 
  dplyr::select(-Physical_pos) %>% 
  dplyr::select(Chromosome, 
                SNP, 
                Genetic_pos, 
                Position)

  

map = read_tsv('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.map', 
               col_names = c('Chromosome', 
                             'SNP', 
                             'Genetic_pos', 
                             'Position'))


inner_join(new_map_clean, 
           map) %>% 
  write_tsv('AC_New_New_map.map')
s