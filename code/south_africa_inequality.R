############################
# Analysis for South Africa
############################

library(fastverse)
set_collapse(nthreads = 4, na.rm = FALSE)

SAIWI <- fread("data/South Africa_estimated_wealth_index.csv") %>% frename(estimated_IWI = IWI)
SARWI <- fread("data/South Africa_RWI.csv") %>% frename(longitude = lon, latitude = lat, rwi = RWI)

source("code/gini_helper.R")
source("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/code/osm_helpers.R")

# South Africa Boundaries
minlon = 10; minlat = -36; maxlon = 40; maxlat = -20 

# Adding WorldPop Gridded Population at 1km Resolution
WPOP <- terra::rast("/Users/sebastiankrantz/Documents/Data/WorldPop/africa_pop_2020_1km.tif")
SAWPOP <- as.data.frame(WPOP, xy = TRUE) %>% qDT() %>% 
  frename(x = lon, y = lat, africa_pop_2020_1km = pop) %>% 
  fsubset(between(lon, minlon, maxlon) & between(lat, minlat, maxlat))

# Adding VIIRS Nightlights, 500m resolution, aggregate by factor 2 to 1km resolution
# NL21MED <- terra::rast("/Users/sebastiankrantz/Documents/Data/VIIRS-DNB/VNL_v21_npp_2021_global_vcmslcfg_c202205302300.median_masked.dat.tif")
# NL21MED %>% terra::crop(terra::ext(c(minlon, maxlon, minlat, maxlat)), filename = "data/SA_VNL_v21_npp_2021_global_vcmslcfg_median_masked.dat.tif") 
SANL21 <- terra::rast("data/SA_VNL_v21_npp_2021_global_vcmslcfg_median_masked.dat.tif") %>% 
          terra::aggregate(fact = 2) %>% as.data.frame(xy = TRUE) %>% set_names(.c(lon, lat, NL21)) %>% 
          fsubset(between(lon, minlon, maxlon) & between(lat, minlat, maxlat))


# Rounding and merging with population
SARWI %$% round_to_kms_exact(lon, lat, 2.1) %>% GRPN(FALSE) %>% descr()
SAIWI %$% round_to_kms_exact(lon, lat, 1.2) %>% GRPN(FALSE) %>% descr()
SANL21 %$% round_to_kms_fast(lon, lat, 0.7) %>% GRPN(FALSE) %>% descr()
SAWPOP %$% round_to_kms_fast(lon, lat, 0.7) %>% GRPN(FALSE) %>% descr()

SARWI_round <- SARWI %>% ftransform(round_to_kms_exact(lon, lat, 2.1)) %>% 
  fgroup_by(lon, lat, sort = FALSE) %>% fselect(RWI) %>% fmean() %>% 
  merge(SAWPOP %>% ftransform(round_to_kms_exact(lon, lat, 2.1)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat))

SAIWI_round <- SAIWI %>% ftransform(round_to_kms_exact(lon, lat, 1.2)) %>% 
  fgroup_by(lon, lat, sort = FALSE) %>% fselect(IWI) %>% fmean() %>% 
  merge(SAWPOP %>% ftransform(round_to_kms_exact(lon, lat, 1.2)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat))

SANL21_round <- SANL21 %>% ftransform(round_to_kms_exact(lon, lat, 0.7)) %>% 
  fgroup_by(lon, lat, sort = FALSE) %>% fselect(NL21) %>% fmean() %>% 
  merge(SAWPOP %>% ftransform(round_to_kms_exact(lon, lat, 0.7)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat))

rm(SANL21, SAWPOP, SAIWI, SARWI); gc()

# 10km Hexagonal grid 
fastverse_extend(qs, dggridR, sf)
qreadm("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/data/intermediate/africa_grids/africa_10km_hexagonal.qs")

# Adding cell identifiers
settransform(SARWI_round, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)
settransform(SAIWI_round, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)
settransform(SANL21_round, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)

# Optimal Power weight for Nightlights to Predict IWI
SA_NL21_IWI <- SAIWI_round %>% 
  merge(SANL21_round %>% ftransform(round_to_kms_exact(lon, lat, 1.2)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fselect(NL21) %>% fsum(), by = .c(lon, lat))

# Optimal Power weight: 0.2749648
opw <- function(a, x, y) cor(x, y^a)
optimise(opw, c(0, 10), x = SA_NL21_IWI$IWI, y = SA_NL21_IWI$NL21, maximum = TRUE)

# Same at grid-level??
settransform(SANL21_round, NL21_scaled = NL21^0.2749648)


# Computing weighted GINI index on Hexagonal grid

SARWI_ineq_10km_hex <- SARWI_round %>% 
  fsubset(GRPN(cell) > 5L) %>% 
  fmutate(RWI = fmin(RWI, TRA = "-", set = TRUE) %+=% 1) %>% 
  fgroup_by(cell) %>% 
  fsummarise(N = GRPN(),
             pop = fsum(pop),
             mean = fmean(RWI),
             w_mean = fmean(RWI, pop),
             median = fmedian(RWI),
             w_median = fmedian(RWI, pop, ties = "q7"), 
             GINI = gini_wiki(RWI), # Gini Index: Area above Lorentz Curve
             w_GINI = w_gini(RWI, pop),
             TI = fmean(fmean(RWI, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(RWI, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index

SAIWI_ineq_10km_hex <- SAIWI_round %>% 
  fsubset(GRPN(cell) > 5L) %>% 
  fgroup_by(cell) %>% 
  fsummarise(N = GRPN(),
             pop = fsum(pop),
             mean = fmean(IWI),
             w_mean = fmean(IWI, pop),
             median = fmedian(IWI),
             w_median = fmedian(IWI, pop, ties = "q7"), 
             GINI = gini_wiki(IWI), # Gini Index: Area above Lorentz Curve
             w_GINI = w_gini(IWI, pop),
             TI = fmean(fmean(IWI, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(IWI, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index

SANL21_ineq_10km_hex <- SANL21_round %>% 
  fsubset(GRPN(cell) > 5L & fsd(NL21_scaled, cell, TRA = 1) > 0) %>% 
  fgroup_by(cell) %>% 
  fsummarise(N = GRPN(),
             pop = fsum(pop),
             mean = fmean(NL21_scaled),
             w_mean = fmean(NL21_scaled, pop),
             median = fmedian(NL21_scaled),
             w_median = fmedian(NL21_scaled, pop, ties = "q7"), 
             GINI = gini_wiki(NL21_scaled), # Gini Index: Area above Lorentz Curve
             w_GINI = w_gini(NL21_scaled, pop),
             TI = fmean(fmean(NL21_scaled+1, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(NL21_scaled, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index


# Merging
SA_ineq_10km_hex_all <- SAIWI_ineq_10km_hex %>% add_stub("IWI_", cols = -1L) %>% 
  merge(SARWI_ineq_10km_hex %>% add_stub("RWI_", cols = -1L), by = "cell", all = TRUE) %>% 
  merge(SANL21_ineq_10km_hex %>% add_stub("NL21_", cols = -1L), by = "cell", all = TRUE) %>% 
  fsubset(is.finite(IWI_mean) | is.finite(RWI_mean) | is.finite(NL21_mean))

fnobs(SA_ineq_10km_hex_all) # check if all matched, but seems so
SA_ineq_10km_hex_all %>% fselect(IWI_pop, RWI_pop) %>% na_omit() %$% cor(IWI_pop, RWI_pop)

# Checking correlations
SA_ineq_10km_hex_all %$% pwcor(IWI_w_GINI, RWI_w_GINI)

# Adding coordinates
add_vars(SA_ineq_10km_hex_all) <- fsubset(africa_10km_hex_sf, ckmatch(SA_ineq_10km_hex_all$cell, cell), -cell)
SA_ineq_10km_hex_all %<>% st_as_sf()

# Saving as geo pakage
st_write(SA_ineq_10km_hex_all, "data/SA_ineq_10km_hex_all.gpkg")