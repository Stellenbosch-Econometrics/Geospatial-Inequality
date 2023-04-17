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
SA_POP20 <- terra::rast("data/WPOP_SA_1km_UNadj/zaf_ppp_2020_1km_Aggregated_UNadj.tif") %>% 
            as.data.frame(xy = TRUE) %>% set_names(.c(lon, lat, POP20)) %>% qDT()

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


add_vars(UberH7) <- ss(H7_grid_all, match(UberH7$hex7, H7_grid_all$hex7), -1L)
st_write(UberH7, "data/spatial_tax_panel/H7_grid_all.gpkg")
rm(UberH7)
H7_grid_all <- st_read("data/spatial_tax_panel/H7_grid_all.gpkg")

# Exploring the grid

H7_grid_all %>% 
  fselect(IWI, RWI, NL20, POP20, MedianIncome_2020, FTE_2020, gini_2020) %>% 
  pwcor(P = TRUE, N = TRUE) %>% print(digits = 3, sig.level = 0.05, show = "lower.tri")

H7_grid_all %>% 
  fmutate(lat_lon_int = finteraction(as.integer(lat), as.integer(lon))) %>% 
  fcomputev(c(IWI, RWI, NL20, MedianIncome_2020, FTE_2020, gini_2020), 
            fscale, lat_lon_int, na.rm = TRUE) %>% 
  pwcor(P = TRUE, N = TRUE) %>% print(digits = 3, sig.level = 0.05, show = "lower.tri")

H7_grid_all %>% 
  fcomputev(c(IWI, RWI, NL20, MedianIncome_2020, FTE_2020, gini_2020), 
            fscale, CAT_B, na.rm = TRUE) %>% 
  pwcor(P = TRUE, N = TRUE) %>% print(digits = 3, sig.level = 0.05, show = "lower.tri")


# Now also getting municipal geometry
MunGeo <- st_read("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Shapefiles/MDB_Local_Municipal_Boundary_2018")
IWI_mun <- SA_IWI %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% num_vars() %>% aggregate(MunGeo, mean) 
RWI_mun <- SA_RWI %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% num_vars() %>% aggregate(MunGeo, mean) 
NL20_mun <- SA_NL20 %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% num_vars() %>% aggregate(MunGeo, mean) 
POP20_mun <- SA_POP20 %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% num_vars() %>% aggregate(MunGeo, sum) 

all_identical(IWI_mun$geometry, RWI_mun$geometry, NL20_mun$geometry, POP20_mun$geometry, MunGeo$geometry)

MunGeo_data <- MunGeo %>% 
  add_vars(nv(unclass(IWI_mun)), nv(unclass(RWI_mun)), nv(unclass(NL20_mun)), nv(unclass(POP20_mun)), 
           sbt(MUN_wide, ckmatch(MunGeo$CAT_B, CAT_B), -1L))

st_write(MunGeo_data, "data/spatial_tax_panel/MunGeo_data.gpkg")
MunGeo_data <- st_read("data/spatial_tax_panel/MunGeo_data.gpkg")

# Plotting 2020 STP3 Gini
range_title <- function(title, x, digits = 2) {
  paste(title, do.call(sprintf, as.list(c("[%s, %s]", round(.range(x), digits)))))
}

plot(MunGeo_data["gini_2020"] %>% mtt(gini_2020 = gini_2020*100), key.pos = 4, border = NA, 
     main = range_title("STP3: GINI Index 2020", MunGeo_data$gini_2020*100))
dev.copy(pdf, "figures/STP3_municipal_gini_2020.pdf", width = 8.27, height = 5.83)
dev.off()

# Plot STP3 GINI for all years
oldpar <- par(mfrow = c(3, 3), oma = c(0,0,0,0), mar = c(0,0,0,0))
for (y in 2014:2022) {
  v = sprintf("gini_%i", y)
  plot(MunGeo_data[v], main = range_title(sprintf("GINI Index %i", y), MunGeo_data[[v]]), 
       border = NA, key.pos = NULL, reset = FALSE, range = c(0.2, 1))
}
par(oldpar)

dev.copy(png, "figures/STP3_municipal_gini.png", width = 11.69*200, height = 8.27*200)
dev.off()


# Correaltion of Nightlights and Other measures 
MunGeo_data %$% cbind(IWI, RWI, log10_NL20 = log10(NL20)) %>% pwcor()

MunGeo_data %>% qDT() %>% 
  fcompute(Log10_NL20 = log10(NL20),
           MedInc_2020 = MedianIncome_2020,
           keep = .c(IWI, RWI, POP20, FTE_2020)) %>% 
  colorder(IWI, RWI, Log10_NL20, POP20, FTE_2020, MedInc_2020) %>% 
  pwcor(P = TRUE) %>% 
  print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE) %>%
  xtable::xtable() %>% print(booktabs = TRUE)

# MunGeo_data %>% fselect(MedianIncome_2020, GINI_2020 = gini_2020, RWI, IWI, NL20, FTE_2020) %>% plot()
oldpar <- par(mfrow = c(2, 3), oma = c(0,0,0,0), mar = c(0,0,0,0))
plot(MunGeo_data["IWI"], main = range_title("International Wealth Index", MunGeo_data$IWI), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["RWI"], main = range_title("Relative Wealth Index", MunGeo_data$RWI), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["NL20"] %>% fcompute(Log10_NL20 = log10(NL20)), main = range_title("Log10 VIIRS-DNB 2020", log10(MunGeo_data$NL20)), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["POP20"], main = range_title("WorldPop Population 2020", MunGeo_data$POP20), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["FTE_2020"], main = range_title("STP3: Full-Time-Eq. Employees 2020", MunGeo_data$FTE_2020), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["MedianIncome_2020"], main = range_title("STP3: Median Income 2020", MunGeo_data$MedianIncome_2020), border = NA, key.pos = NULL, reset = FALSE)
par(oldpar)

dev.copy(pdf, "figures/STP3_municipal_joint.pdf", width = 11.69, height = 8.27)
dev.off()

#
### Matching municipalities and Uber Hexagons ---------
#
set_collapse(na.rm = TRUE)
MunGeo_data <- st_read("data/spatial_tax_panel/MunGeo_data.gpkg")
H7_grid_all <- st_read("data/spatial_tax_panel/H7_grid_all.gpkg")

ind = st_contains(MunGeo_data$geom, st_centroid(H7_grid_all$geom))
H7_grid_all$CAT_B <- NA_character_
for(i in seq_along(ind)) setv(H7_grid_all$CAT_B, ind[[i]], MunGeo_data$CAT_B[i], vind1 = TRUE)
# Now regenerate H7_grid_all
# Saving
st_write(H7_grid_all, "data/spatial_tax_panel/H7_grid_all.gpkg")

# Computing GINI estimates
source("code/gini_helper.R")


H7_grid_all %>% fselect(IWI, RWI, NL20, POP20) %>% descr()

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
    fsummarise(NL20_W_GINI = w_gini(NL20^k, POP20, na.rm = TRUE)) %>% # gini_wiki(na_rm(NL20^k))
    merge(Mun_gini_19_21_avg, by = "CAT_B") %$%
    fmean(abs(NL20_W_GINI - gini_avg))
}
NL20_opt <- optimise(NL20_objective, c(0.0001, 100))
# Unweighted: k = 1.566813, Obj = 0.1866354

IWI_objective <- function(k) {
  H7_grid_all_grouped %>% 
    fsummarise(IWI_W_GINI = w_gini(IWI^k, POP20, na.rm = TRUE)) %>% # gini_wiki(na_rm(IWI^k))
    merge(Mun_gini_19_21_avg, by = "CAT_B") %$%
    fmean(abs(IWI_W_GINI - gini_avg))
}
IWI_opt <- optimise(IWI_objective, c(0.0001, 100))
# Unweighted: k = 6.537822, Obj = 0.1288476 (better!!)

RWI_objective <- function(k) {
  H7_grid_all_grouped %>% 
    fsummarise(RWI_W_GINI = w_gini(RWI^k, POP20, na.rm = TRUE)) %>% # gini_wiki(na_rm(RWI^k))
    merge(Mun_gini_19_21_avg, by = "CAT_B") %$%
    fmean(abs(RWI_W_GINI - gini_avg))
}
RWI_opt <- optimise(RWI_objective, c(0.0001, 100))
# Unweighted: k = 4.070419, Obj = 0.08430859 (better!!)

Mun_GINI_raw <- H7_grid_all_grouped %>% 
  fsummarise(NL20_GINI = gini_wiki(na_rm(NL20)), 
             NL20_W_GINI = w_gini(NL20, POP20, na.rm = TRUE),
             IWI_GINI = gini_wiki(na_rm(IWI)),
             IWI_W_GINI = w_gini(IWI, POP20, na.rm = TRUE),
             RWI_GINI = gini_wiki(na_rm(RWI)),
             RWI_W_GINI = w_gini(RWI, POP20, na.rm = TRUE)) %>% 
  merge(Mun_gini_19_21_avg, by = "CAT_B")

Mun_GINI_calib <- H7_grid_all_grouped %>% 
  fsummarise(NL20_GINI = gini_wiki(na_rm(NL20^1.566813)), 
             NL20_W_GINI = w_gini(NL20^0.1791315, POP20, na.rm = TRUE),
             IWI_GINI = gini_wiki(na_rm(IWI^6.537822)),
             IWI_W_GINI = w_gini(IWI^14.71927, POP20, na.rm = TRUE),
             RWI_GINI = gini_wiki(na_rm(RWI^4.070419)),
             RWI_W_GINI = w_gini(RWI^7.044468, POP20, na.rm = TRUE)) %>% 
  merge(Mun_gini_19_21_avg, by = "CAT_B")


Mun_GINI_raw %>% num_vars() %>% pwcor()
Mun_GINI_calib %>% num_vars() %>% pwcor()
             


 
