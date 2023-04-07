###########################
# Explore Spatial Tax Panel
###########################

library(fastverse)
set_collapse(nthreads = 4, na.rm = FALSE)

STP3_files <- list.files("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data", full.names = TRUE)
names(STP3_files) <- flast(strsplit(STP3_files, "/", fixed = TRUE)) %>% sub(pat = ".csv", rep = "")
STP3_Municipal <- STP3_files[names(STP3_files) %ilike% "Gini|MedianIncome" & startsWith(names(STP3_files), "Municipal")] %>% lapply(fread)
STP3_H7 <- STP3_files[names(STP3_files) %ilike% "Gini|MedianIncome" & startsWith(names(STP3_files), "hex7")] %>% lapply(fread)

# FTE = The number of full time equivalent (FTE) employees.
# gini = The gini coefficient calculated based on the weighted incomes of all employees within an aggregation.

# Spatial Tax Panel: Municipal
MUN <- fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/Municipal_MedianIncome.csv") %>% fselect(-FTE) %>% 
       merge(fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/Municipal_Gini.csv"), by = .c(CAT_B, TaxYear)) %>% 
       fmutate(FTE = as.double(FTE))

descr(MUN)
hist(MUN$gini, breaks = 100)

MUN_wide <- MUN %>% dcast(CAT_B ~ TaxYear, value.var = .c(MedianIncome, FTE, gini)) %T>% 
  fwrite("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/Municipal_MedianIncome_Gini_wide.csv")

# Spatial Tax Panel: Hexagons
H7 <- fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/hex7_MedianIncome.csv") %>% fselect(-FTE, -CAT_B) %>% 
      merge(fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/hex7_Gini.csv"), by = .c(hex7, TaxYear)) %>% 
      fmutate(FTE = as.double(FTE))

descr(H7)
hist(H7$gini, breaks = 100)

H7_wide <- H7 %>% dcast(hex7 ~ TaxYear, value.var = .c(MedianIncome, FTE, gini)) %T>% 
  fwrite("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/hex7_MedianIncome_Gini_wide.csv")

# Getting H7 Geometry
fastverse_extend(sf)
source("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/code/osm_helpers.R")

UberH7 <- st_read("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Shapefiles/UberH3_7")
UberH7_centroids <- UberH7 %>% st_centroid() %>% qDT() %>% 
  ftransform(st_coordinates(geometry) %>% mctl(TRUE) %>% set_names(c("lon", "lat"))) %>% 
  ftransform(geometry = NULL)
# Spatial resolution: 1.8km  
UberH7_centroids %$% round_to_kms_exact(lon, lat, 1.8) %>% any_duplicated()

# Adding Wealth Estimates
SA_IWI <- fread("data/South Africa_estimated_wealth_index.csv") %>% frename(estimated_IWI = IWI) %>% fmutate(img_embedding = NULL)
SA_RWI <- fread("data/South Africa_RWI.csv") %>% frename(longitude = lon, latitude = lat, rwi = RWI)
SA_NL20 <- terra::rast("data/south_africa_viirs_dnb_nightlights_v1_vcmslcfg_annual_median_composite/2020.tif") %>% 
           as.data.frame(xy = TRUE) %>% set_names(.c(lon, lat, NL20)) %>% fsubset(NL20 %!=% -9999) %>% qDT()

# Merging
H7_grid_all <- UberH7_centroids %>% ftransform(round_to_kms_exact(lon, lat, 1.8)) %>%
  merge(SA_IWI %>% ftransform(round_to_kms_exact(lon, lat, 1.8)) %>% 
        fgroup_by(lon, lat) %>% num_vars() %>% fmean(), by = .c(lon, lat), all.x = TRUE) %>% 
  merge(SA_RWI %>% ftransform(round_to_kms_exact(lon, lat, 1.8)) %>% 
        fgroup_by(lon, lat) %>% num_vars() %>% fmean(), by = .c(lon, lat), all.x = TRUE) %>% 
  merge(SA_NL20 %>% ftransform(round_to_kms_exact(lon, lat, 1.8)) %>% 
        fgroup_by(lon, lat) %>% num_vars() %>% fmean(), by = .c(lon, lat), all.x = TRUE) %>% 
  merge(H7_wide, by = "hex7", all.x = TRUE) %>% 
  na_omit(cols = -(1:3), prop = 1) %>% roworder(hex7)

# Exploring the grid
H7_grid_all %>% fselect(IWI, RWI, NL20, MedianIncome_2020, FTE_2020, gini_2020) %>% pwcor()

