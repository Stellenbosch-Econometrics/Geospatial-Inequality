###########################################
# Population Weighted Geospatial Inequality
###########################################

library(fastverse)
fastverse_extend(qs) 
set_collapse(nthreads = 4, na.rm = FALSE)

source("code/gini_helper.R")
source("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/code/osm_helpers.R")

# Import data and cropping to Africa
minlon = -27; minlat = -36; maxlon = 59; maxlat = 38 # Africa Boundaries

RWI <- qread("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/data/intermediate/other/RWI.qs") %>% 
  fsubset(between(lon, minlon, maxlon) & between(lat, minlat, maxlat))
descr(RWI)
IWI <- qread("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/data/intermediate/other/IWI.qs") %>% 
  fsubset(between(lon, minlon, maxlon) & between(lat, minlat, maxlat))
descr(IWI)
# NL21 <- qread("/Users/sebastiankrantz/Documents/Data/VIIRS-DNB/viirs_dnb_2021_median_masked.qs") %>% 
#   fsubset(between(lon, minlon, maxlon) & between(lat, minlat, maxlat)) %>% qDT()
# qsu(NL21, stable.algo = FALSE)

# Adding WorldPop Gridded Population at 1km Resolution
fastverse_extend(terra)
WPOP <- rast("/Users/sebastiankrantz/Documents/Data/WorldPop/africa_pop_2020_1km.tif")
# plot(WPOP)
WPOP <- as.data.frame(WPOP, xy = TRUE) %>% qDT()
setrename(WPOP, x = lon, y = lat, africa_pop_2020_1km = pop)

# Matching Population to RWI and IWI ---------------------------------------------

RWI %$% round_to_kms_exact(lon, lat, 2.4) %>% GRPN(FALSE) %>% descr()
IWI %$% round_to_kms_exact(lon, lat, 1.7) %>% GRPN(FALSE) %>% descr()

RWI_round <- RWI %>% ftransform(round_to_kms_exact(lon, lat, 2.4)) %>% 
      fgroup_by(lon, lat, sort = FALSE) %>% fselect(RWI) %>% fmean() %>% 
      merge(WPOP %>% ftransform(round_to_kms_exact(lon, lat, 2.4)) %>% 
              fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat))

IWI_round <- IWI %>% ftransform(round_to_kms_exact(lon, lat, 1.7)) %>% 
  fgroup_by(lon, lat, sort = FALSE) %>% fselect(IWI) %>% fmean() %>% 
  merge(WPOP %>% ftransform(round_to_kms_exact(lon, lat, 1.7)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat))


# 10km Hexagonal grid 
fastverse_extend(dggridR, sf)
qreadm("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/data/intermediate/africa_grids/africa_10km_hexagonal.qs")

settransform(RWI_round, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)
settransform(IWI_round, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)

# Computing weighted GINI index on Hexagonal grid

RWI_WGINI_10km_hex <- RWI_round %>% 
  fsubset(GRPN(cell) > 5L) %>% 
  fmutate(RWI = fmin(RWI, TRA = "-", set = TRUE) %+=% 1) %>% 
  fgroup_by(cell) %>% 
  fsummarise(N = GRPN(),
             pop = fsum(pop),
             mean = fmean(RWI, pop),
             median = fmedian(RWI, pop, ties = "q7"), 
             GINI = w_gini(RWI, pop))
             
IWI_WGINI_10km_hex <- IWI_round %>% 
  fsubset(GRPN(cell) > 5L) %>% 
  fgroup_by(cell) %>% 
  fsummarise(N = GRPN(),
             pop = fsum(pop),
             mean = fmean(IWI, pop),
             median = fmedian(IWI, pop, ties = "q7"), 
             GINI = w_gini(IWI, pop))

# Merging
WGINI_10km_hex_all <- IWI_WGINI_10km_hex %>% add_stub("IWI_", cols = -1L) %>% 
  merge(RWI_WGINI_10km_hex %>% add_stub("RWI_", cols = -1L), by = "cell", all = TRUE) %>% 
  fsubset(is.finite(IWI_mean) | is.finite(RWI_mean))

fnobs(WGINI_10km_hex_all) # check if all matched, but seems so
WGINI_10km_hex_all %>% fselect(IWI_pop, RWI_pop) %>% na_omit() %$% cor(IWI_pop, RWI_pop)

# Checking correlations
WGINI_10km_hex_all %$% pwcor(IWI_GINI, RWI_GINI)

# Adding coordinates
add_vars(WGINI_10km_hex_all) <- fsubset(africa_10km_hex_sf, ckmatch(WGINI_10km_hex_all$cell, cell), -cell)
WGINI_10km_hex_all %<>% st_as_sf()

# Saving as geo pakage
st_write(WGINI_10km_hex_all, "data/WGINI_10km_hex_all.gpkg")
