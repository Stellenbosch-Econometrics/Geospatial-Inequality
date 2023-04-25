###########################
# Spatial Tax Panel v3.7
###########################

library(fastverse)
fastverse_extend(sf, install = TRUE)
source("code/spatial_helpers.R")
source("code/gini.R")
set_collapse(nthreads = 4)

#
# Loading all data of interest ----------------------------------------------------------------
# 

STP3_files <- list.files("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data", full.names = TRUE)
names(STP3_files) <- flast(strsplit(STP3_files, "/", fixed = TRUE)) %>% sub(pat = ".csv", rep = "")
STP3_Municipal <- STP3_files[names(STP3_files) %ilike% "Gini$|MedianIncome$" & startsWith(names(STP3_files), "Municipal")] %>% lapply(fread)
STP3_H7 <- STP3_files[names(STP3_files) %ilike% "Gini$|MedianIncome$" & startsWith(names(STP3_files), "hex7")] %>% lapply(fread)
# FTE = The number of full time equivalent (FTE) employees.
# gini = The gini coefficient calculated based on the weighted incomes of all employees within an aggregation.

# Municipal geometry (need to unzip first)
MunGeo_file <- "data/spatial_tax_panel/Spatial_Tax_Panel_v3/Shapefiles/MDB_Local_Municipal_Boundary_2018"
unzip(paste0(MunGeo_file, ".zip"), exdir = MunGeo_file, junkpaths = TRUE)
MunGeo <- st_read(MunGeo_file)

# Uber H3.7 Geometry
UberH7_file <- "data/spatial_tax_panel/Spatial_Tax_Panel_v3/Shapefiles/UberH3_7"
unzip(paste0(UberH7_file, ".zip"), exdir = UberH7_file, junkpaths = TRUE)
UberH7 <- st_read(UberH7_file)

# Load Wealth Estimates
SA_IWI <- fread("data/SA_IWI.csv") %>% frename(estimated_IWI = IWI)
SA_RWI <- fread("data/SA_RWI.csv") %>% frename(longitude = lon, latitude = lat, rwi = RWI)
SA_NL20 <- terra::rast("data/south_africa_viirs_dnb_nightlights_v1_vcmslcfg_annual_median_composite/2020.tif") %>% 
  as.data.frame(xy = TRUE) %>% set_names(.c(lon, lat, NL20)) %>% fsubset(NL20 %!=% -9999) %>% qDT()
SA_POP20 <- terra::rast("data/WPOP_SA_1km_UNadj/zaf_ppp_2020_1km_Aggregated_UNadj.tif") %>% 
  as.data.frame(xy = TRUE) %>% set_names(.c(lon, lat, POP20)) %>% qDT()


#
# Municipalities --------------------------------------------------------------------------
# 

MUN <- STP3_Municipal$Municipal_MedianIncome %>% fselect(-FTE) %>% 
       merge(STP3_Municipal$Municipal_Gini, by = .c(CAT_B, TaxYear)) %>% 
       fmutate(FTE = as.double(FTE))

descr(MUN)
hist(MUN$gini, breaks = 100)

# Reshape wide
MUN_wide <- MUN %>% dcast(CAT_B ~ TaxYear, value.var = .c(MedianIncome, FTE, gini))

# Aggregating data over geometry
IWI_mun <- SA_IWI %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% num_vars() %>% aggregate(MunGeo, mean) 
RWI_mun <- SA_RWI %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% num_vars() %>% aggregate(MunGeo, mean) 
NL20_mun <- SA_NL20 %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% num_vars() %>% aggregate(MunGeo, mean) 
POP20_mun <- SA_POP20 %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% num_vars() %>% aggregate(MunGeo, sum) 

# Check
all_identical(IWI_mun$geometry, RWI_mun$geometry, NL20_mun$geometry, POP20_mun$geometry, MunGeo$geometry)

# Adding aggregated data
MunGeo_data <- MunGeo %>% 
  add_vars(nv(unclass(IWI_mun)), nv(unclass(RWI_mun)), nv(unclass(NL20_mun)), nv(unclass(POP20_mun)), 
           sbt(MUN_wide, ckmatch(MunGeo$CAT_B, CAT_B), -1L))


#
# Uber H3.7 Hexagons -------------------------------------------------------------------------
# 

H7 <- STP3_H7$hex7_MedianIncome %>% fselect(-FTE, -CAT_B) %>% 
      merge(STP3_H7$hex7_Gini, by = .c(hex7, TaxYear)) %>% 
      fmutate(FTE = as.double(FTE))

descr(H7)
hist(H7$gini, breaks = 100)

H7_wide <- H7 %>% dcast(hex7 ~ TaxYear, value.var = .c(MedianIncome, FTE, gini)) 

UberH7_centroids <- UberH7 %>% st_centroid() %>% qDT() %>% 
  ftransform(st_coordinates(geometry) %>% mctl(TRUE) %>% set_names(c("lon", "lat"))) %>% 
  ftransform(geometry = NULL)

# Spatial resolution: 1.8km  
UberH7_centroids %$% round_to_kms_exact(lon, lat, 1.8) %>% any_duplicated()

# Merging
H7_grid_all <- UberH7_centroids %>% ftransform(round_to_kms_exact(lon, lat, 1.8)) %>%
  merge(SA_IWI %>% ftransform(round_to_kms_exact(lon, lat, 1.8)) %>% 
          fgroup_by(lon, lat) %>% num_vars() %>% fmean(), by = .c(lon, lat), all.x = TRUE) %>% 
  merge(SA_RWI %>% ftransform(round_to_kms_exact(lon, lat, 1.8)) %>% 
          fgroup_by(lon, lat) %>% num_vars() %>% fmean(), by = .c(lon, lat), all.x = TRUE) %>% 
  merge(SA_NL20 %>% ftransform(round_to_kms_exact(lon, lat, 1.8)) %>% 
          fgroup_by(lon, lat) %>% num_vars() %>% fmean(), by = .c(lon, lat), all.x = TRUE) %>% 
  merge(SA_POP20 %>% ftransform(round_to_kms_exact(lon, lat, 1.8)) %>% 
          fgroup_by(lon, lat) %>% num_vars() %>% fsum(), by = .c(lon, lat), all.x = TRUE) %>% 
  merge(H7_wide, by = "hex7", all.x = TRUE) %>% 
  na_omit(cols = -(1:3), prop = 1) %>% roworder(hex7)

H7_grid_all %>% fselect(IWI, RWI, NL20, POP20) %>% descr()

# Combining and renaming
add_vars(UberH7) <- ss(H7_grid_all, match(UberH7$hex7, H7_grid_all$hex7), -1L)
H7_grid_all <- UberH7
rm(UberH7)

# Matching hexagons to municipalities
ind <- st_contains(MunGeo_data$geom, st_centroid(H7_grid_all$geom))
H7_grid_all$CAT_B <- NA_character_
for(i in seq_along(ind)) setv(H7_grid_all$CAT_B, ind[[i]], MunGeo_data$CAT_B[i], vind1 = TRUE)
descr(H7_grid_all$CAT_B) # Note: Some NA's
setcolorder(H7_grid_all, .c(hex7, CAT_B))

# Saving
st_write(H7_grid_all, "results/STP3_Hex7_grid.gpkg")



#
# Computing GINI estimates: calibration to municipal GINIs ----------------------------------------------------
# 

# Grouping Hexagons by Municipality: use the hexagons as baseline for municipal estimates
H7_grid_all_grouped <- H7_grid_all %>% qDT() %>% 
  fsubset(POP20 > 0) %>% 
  fmutate(RWI = RWI - fmin(RWI) + 0.1) %>%
  # fselect(IWI, RWI, NL20, POP20) %>% descr()
  fgroup_by(CAT_B) 


# Target
Mun_gini_19_21_avg <- MunGeo_data %>% qDT() %>% 
  fcompute(gini_avg = pmean(gini_2019, gini_2020, gini_2021, na.rm = TRUE), 
           keep = "CAT_B")

# Now the Optimization
NL20_objective <- function(k) {
  H7_grid_all_grouped %>% 
    fsummarise(NL20_W_GINI = w_gini(NL20^k, POP20, na.rm = TRUE)) %>% # Unweighted: gini_wiki(na_rm(NL20^k))
    merge(Mun_gini_19_21_avg, by = "CAT_B") %$%
    fmean(abs(NL20_W_GINI - gini_avg))
}
NL20_opt <- optimise(NL20_objective, c(0.0001, 100))
# Unweighted: k = 1.566813, Obj = 0.1866354

IWI_objective <- function(k) {
  H7_grid_all_grouped %>% 
    fsummarise(IWI_W_GINI = w_gini(IWI^k, POP20, na.rm = TRUE)) %>% # Unweighted: gini_wiki(na_rm(IWI^k))
    merge(Mun_gini_19_21_avg, by = "CAT_B") %$%
    fmean(abs(IWI_W_GINI - gini_avg))
}
IWI_opt <- optimise(IWI_objective, c(0.0001, 100))
# Unweighted: k = 6.537822, Obj = 0.1288476 (better!!)

RWI_objective <- function(k) {
  H7_grid_all_grouped %>% 
    fsummarise(RWI_W_GINI = w_gini(RWI^k, POP20, na.rm = TRUE)) %>% # Unweighted: gini_wiki(na_rm(RWI^k))
    merge(Mun_gini_19_21_avg, by = "CAT_B") %$%
    fmean(abs(RWI_W_GINI - gini_avg))
}
RWI_opt <- optimise(RWI_objective, c(0.0001, 100))
# Unweighted: k = 4.070419, Obj = 0.08430859 (better!!)

# Raw (uncalibrated) Municipal GINI estmates
Mun_GINI_raw <- H7_grid_all_grouped %>% 
  fsummarise(NL20_GINI = gini_wiki(na_rm(NL20)), 
             NL20_W_GINI = w_gini(NL20, POP20, na.rm = TRUE),
             IWI_GINI = gini_wiki(na_rm(IWI)),
             IWI_W_GINI = w_gini(IWI, POP20, na.rm = TRUE),
             RWI_GINI = gini_wiki(na_rm(RWI)),
             RWI_W_GINI = w_gini(RWI, POP20, na.rm = TRUE)) %>% 
  merge(Mun_gini_19_21_avg, by = "CAT_B")

Mun_GINI_raw %>% num_vars() %>% pwcor()

# Calibrated GINI estimates
Mun_GINI_calib <- H7_grid_all_grouped %>% 
  fsummarise(NL20_GINI = gini_wiki(na_rm(NL20^1.566813)), 
             NL20_W_GINI = w_gini(NL20^1.655078, POP20, na.rm = TRUE),
             IWI_GINI = gini_wiki(na_rm(IWI^6.537822)),
             IWI_W_GINI = w_gini(IWI^14.71927, POP20, na.rm = TRUE),
             RWI_GINI = gini_wiki(na_rm(RWI^4.070419)),
             RWI_W_GINI = w_gini(RWI^7.044468, POP20, na.rm = TRUE)) %>% 
  merge(Mun_gini_19_21_avg, by = "CAT_B")

Mun_GINI_calib %>% num_vars() %>% pwcor()

# Adding estimates to the data
tfm(MunGeo_data) <- Mun_GINI_raw %>% fsubset(match(MunGeo_data$CAT_B, CAT_B), -CAT_B) %c*% 100
tfm(MunGeo_data) <- Mun_GINI_calib %>% fsubset(match(MunGeo_data$CAT_B, CAT_B), -CAT_B, -gini_avg) %>% add_stub("_calib", FALSE) %c*% 100

# Saving results
st_write(MunGeo_data, "results/STP3_MunGeo_data.gpkg")
