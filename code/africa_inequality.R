########################
# Geospatial Inequality
########################

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
NL21 <- qread("/Users/sebastiankrantz/Documents/Data/VIIRS-DNB/viirs_dnb_2021_median_masked.qs") %>% 
       fsubset(between(lon, minlon, maxlon) & between(lat, minlat, maxlat)) %>% qDT()
qsu(NL21, stable.algo = FALSE)


# 20km Rectilinear Grid ---------------------------------------------------------

# Computing inequality measures (rectilinear grid)
RWI_ineq_20km <- RWI %>% 
  # fsubset(ISO3 %==% "ZAF") %>% 
  ftransform(round_to_kms_exact(lon, lat, 20)) %>% 
  # fcount(lon, lat) %$% descr(N)
  fsubset(GRPN(list(lon, lat)) > 5L) %>% 
  fmutate(RWI = fmin(RWI, TRA = "-", set = TRUE) %+=% 1) %>% 
  fgroup_by(lon, lat) %>% 
  fsummarise(ISO3 = fmode(ISO3), 
             N = GRPN(),
             mean = fmean(RWI),
             median = fmedian(RWI), 
             GINI = gini_wiki(RWI), # Gini Index: Area above Lorentz Curve
             TI = fmean(fmean(RWI, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(RWI, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index

IWI_ineq_20km <- IWI %>% 
  ftransform(round_to_kms_exact(lon, lat, 20)) %>% 
  fsubset(GRPN(list(lon, lat)) > 5L) %>% 
  fgroup_by(lon, lat) %>% 
  fsummarise(country_name = fmode(country_name), 
             N = GRPN(),
             mean = fmean(IWI),
             median = fmedian(IWI), 
             GINI = gini_wiki(IWI), # Gini Index: Area above Lorentz Curve
             TI = fmean(fmean(IWI, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(IWI, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index

NL21_ineq_20km <- NL21 %>% 
  ftransform(round_to_kms_exact(lon, lat, 20)) %>% 
  # fcount(lon, lat) %$% descr(N)
  fsubset(GRPN(list(lon, lat)) > 5L) %>% 
  fmutate(nl_21_med_mask = fmin(nl_21_med_mask, TRA = "-", set = TRUE) %+=% 1) %>% 
  fgroup_by(lon, lat) %>% 
  fsummarise(N = GRPN(),
             mean = fmean(nl_21_med_mask),
             median = fmedian(nl_21_med_mask), 
             GINI = gini_wiki(nl_21_med_mask), # Gini Index: Area above Lorentz Curve
             TI = fmean(fmean(nl_21_med_mask, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(nl_21_med_mask, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index

# Merging
ineq_20km_all <- NL21_ineq_20km %>% add_stub("NL21_", cols = -(1:2)) %>% 
  merge(IWI_ineq_20km %>% add_stub("IWI_", cols = -(1:2)), by = c("lon", "lat"), all.x = TRUE) %>% 
  merge(RWI_ineq_20km %>% add_stub("RWI_", cols = -(1:2)), by = c("lon", "lat"), all.x = TRUE) %>% 
  roworder(lon, lat)

fnobs(ineq_20km_all) # check if all matched, but seems so

# Saving 
fwrite(ineq_20km_all, "data/ineq_20km_all.csv")

ineq_20km_all %>% 
  gvr("GINI|NHI|TI") %>% 
  pwcor()


# 10km Hexagonal grid ------------------------------------------------------------

fastverse_extend(dggridR, sf)
qreadm("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/data/intermediate/africa_grids/africa_10km_hexagonal.qs")

settransform(RWI, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)
settransform(IWI, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)
settransform(NL21, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)

# Computing inequality measures (hexagonal grid)
RWI_ineq_10km_hex <- RWI %>% 
  fsubset(GRPN(cell) > 5L) %>% 
  fmutate(RWI = fmin(RWI, TRA = "-", set = TRUE) %+=% 1) %>% 
  fgroup_by(cell) %>% 
  fsummarise(ISO3 = fmode(ISO3), 
             N = GRPN(),
             mean = fmean(RWI),
             median = fmedian(RWI), 
             GINI = gini_wiki(RWI), # Gini Index: Area above Lorentz Curve
             TI = fmean(fmean(RWI, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(RWI, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index

IWI_ineq_10km_hex <- IWI %>% 
  fsubset(GRPN(cell) > 5L) %>% 
  fgroup_by(cell) %>% 
  fsummarise(country_name = fmode(country_name), 
             N = GRPN(),
             mean = fmean(IWI),
             median = fmedian(IWI), 
             GINI = gini_wiki(IWI), # Gini Index: Area above Lorentz Curve
             TI = fmean(fmean(IWI, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(IWI, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index

NL21_ineq_10km_hex <- NL21 %>% 
  fsubset(GRPN(cell) > 5L) %>% 
  fmutate(nl_21_med_mask = fmin(nl_21_med_mask, TRA = "-", set = TRUE) %+=% 1) %>% 
  fgroup_by(cell) %>% 
  fsummarise(N = GRPN(),
             mean = fmean(nl_21_med_mask),
             median = fmedian(nl_21_med_mask), 
             GINI = gini_wiki(nl_21_med_mask), # Gini Index: Area above Lorentz Curve
             TI = fmean(fmean(nl_21_med_mask, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(nl_21_med_mask, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index


# Merging
ineq_10km_hex_all <- NL21_ineq_10km_hex %>% add_stub("NL21_", cols = -1L) %>% 
  merge(IWI_ineq_10km_hex %>% add_stub("IWI_", cols = -1L), by = "cell", all.x = TRUE) %>% 
  merge(RWI_ineq_10km_hex %>% add_stub("RWI_", cols = -1L), by = "cell", all.x = TRUE) %>% 
  fsubset(is.finite(IWI_mean) | is.finite(RWI_mean) | NL21_GINI > 0)

fnobs(ineq_10km_hex_all) # check if all matched, but seems so

# Adding coordinates
add_vars(ineq_10km_hex_all) <- fsubset(africa_10km_hex_sf, ckmatch(ineq_10km_hex_all$cell, cell), -cell)
ineq_10km_hex_all %<>% st_as_sf()

# Adding WorldPop Gridded Population at 1km Resolution
fastverse_extend(terra)
WPOP <- rast("/Users/sebastiankrantz/Documents/Data/WorldPop/africa_pop_2020_1km.tif")
plot(WPOP)
WPOP <- as.data.frame(WPOP, xy = TRUE) %>% qDT()
settransform(WPOP, cell = dgGEO_to_SEQNUM(world_10km_hex, x, y)$seqnum)
WPOP_10km_hex <- WPOP %>% collap(africa_pop_2020_1km ~ cell, fsum, sort = FALSE)

ineq_10km_hex_all$POP20 <- WPOP_10km_hex %>% 
    fsubset(cell %in% ineq_10km_hex_all$cell) %>% 
    fsubset(match(ineq_10km_hex_all$cell, cell), africa_pop_2020_1km) %>% .subset2(1L)

descr(ineq_10km_hex_all$POP20)
descr(log(ineq_10km_hex_all$POP20))

# Creating Gini-multiplied measures
tfm(ineq_10km_hex_all) <- ineq_10km_hex_all %>% qDT() %>% 
  gvr("GINI|NHI|TI") %c*% log(ineq_10km_hex_all$POP20) %>% 
  replace_Inf() %>% add_stub("_logPOP20", pre = FALSE)

# Saving as geo pakage
st_write(ineq_10km_hex_all, "data/ineq_10km_hex_all.gpkg")
