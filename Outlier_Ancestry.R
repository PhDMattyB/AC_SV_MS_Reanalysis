##############################
## Outlier ancestry 
##
## Matt Brachmann (PhDMattyB)
##
## 13.04.2023
##
##############################

# Make --keep files for plink -------------------------------

ped_test = read_table2('Charr_Poly_All_Fixed_coords_maf05_geno95_notbed.ped', 
                       col_names = F)

ped_ids = read_table2('Charr_Poly_All_Fixed_coords_maf05_geno95.fam', 
                      col_names = F) %>%
  dplyr::select(X1,
                X2)

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

##
