############################
# Analysis for South Africa
############################

library(fastverse)
set_collapse(nthreads = 4, na.rm = FALSE)

IWI <- fread("data/SA_IWI.csv") %>% frename(estimated_IWI = IWI)
RWI <- fread("data/SA_RWI.csv") %>% frename(longitude = lon, latitude = lat, rwi = RWI)

source("code/gini.R")
source("code/spatial_helpers.R")

# South Africa Boundaries
minlon = 10; minlat = -36; maxlon = 40; maxlat = -20 

#
# Loading and Preparing Data ---------------------------------------------------------------------------
#

# Adding WorldPop Gridded Population at 1km Resolution
WPOP <- terra::rast("data/WPOP_SA_1km_UNadj/zaf_ppp_2020_1km_Aggregated_UNadj.tif") %>% 
     as.data.frame(WPOP, xy = TRUE) %>% qDT() %>% 
     set_names(.c(lon, lat, pop)) %>% 
     fsubset(between(lon, minlon, maxlon) & between(lat, minlat, maxlat))

# Adding VIIRS Nightlights, 463.8m resolution, aggregate by factor 2 to ~1km resolution
NL21 <- terra::rast("data/SA_VNL_v21/SA_VNL_v21_npp_2021_global_vcmslcfg_median_masked.dat.tif") %>%
     terra::aggregate(fact = 2, fun = "mean", na.rm = TRUE) %>% 
     as.data.frame(xy = TRUE) %>% set_names(.c(lon, lat, NL21)) %>% 
     fsubset(between(lon, minlon, maxlon) & between(lat, minlat, maxlat))

# Rounding and merging with population
RWI %$% round_to_kms_exact(lon, lat, 2.1) %>% GRPN(FALSE) %>% descr()
IWI %$% round_to_kms_exact(lon, lat, 1.2) %>% GRPN(FALSE) %>% descr()
NL21 %$% round_to_kms_fast(lon, lat, 0.7) %>% GRPN(FALSE) %>% descr()
WPOP %$% round_to_kms_fast(lon, lat, 0.7) %>% GRPN(FALSE) %>% descr()

RWI_round <- RWI %>% ftransform(round_to_kms_exact(lon, lat, 2.1)) %>% 
  fgroup_by(lon, lat, sort = FALSE) %>% fselect(RWI) %>% fmean() %>% 
  merge(WPOP %>% ftransform(round_to_kms_exact(lon, lat, 2.1)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat)) %>% 
  fmutate(RWI = fmin(RWI, TRA = "-", set = TRUE) %+=% 1)

IWI_round <- IWI %>% ftransform(round_to_kms_exact(lon, lat, 1.2)) %>% 
  fgroup_by(lon, lat, sort = FALSE) %>% fselect(IWI) %>% fmean() %>% 
  merge(WPOP %>% ftransform(round_to_kms_exact(lon, lat, 1.2)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat))

NL21_round <- NL21 %>% ftransform(round_to_kms_exact(lon, lat, 0.7)) %>% 
  fgroup_by(lon, lat, sort = FALSE) %>% fselect(NL21) %>% fmean() %>% 
  merge(WPOP %>% ftransform(round_to_kms_exact(lon, lat, 0.7)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat))

# Removing original files and freeing memory (garbage collection)
rm(NL21, WPOP, IWI, RWI); gc()


#
# Inequality Within 96km2 Hexagonal Grid -----------------------------------------------------------------------------------
#

# Loading / installing spatial libraries
fastverse_extend(dggridR, sf, install = TRUE)

# This constructs a uniform 96km2 hexagonal grid of the Earth's Surface
world_96km2_hex <- dgconstruct(area = 144, metric = TRUE, resround = "nearest")

# Adding grid cell identifiers
settransform(RWI_round, cell = dgGEO_to_SEQNUM(world_96km2_hex, lon, lat)$seqnum)
settransform(IWI_round, cell = dgGEO_to_SEQNUM(world_96km2_hex, lon, lat)$seqnum)
settransform(NL21_round, cell = dgGEO_to_SEQNUM(world_96km2_hex, lon, lat)$seqnum)

# Computing (weighted) GINI, Theil Index, and (normalized) Herfindahl Index on Hexagonal Grid

RWI_ineq_96km2_hex <- RWI_round %>% 
  fsubset(GRPN(cell) > 5L) %>% 
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

IWI_ineq_96km2_hex <- IWI_round %>% 
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

NL21_ineq_96km2_hex <- NL21_round %>% 
  fsubset(GRPN(cell) > 5L & fsd(NL21, cell, TRA = 1) > 0) %>% 
  fgroup_by(cell) %>% 
  fsummarise(N = GRPN(),
             pop = fsum(pop),
             mean = fmean(NL21),
             w_mean = fmean(NL21, pop),
             median = fmedian(NL21),
             w_median = fmedian(NL21, pop, ties = "q7"), 
             GINI = gini_wiki(NL21), # Gini Index: Area above Lorentz Curve
             w_GINI = w_gini(NL21, pop),
             TI = fmean(fmean(NL21+1, TRA = "/") %>% multiply_by(log(.))), # Theil Index (Overall)                                                    
             HI = fsum(fsum(NL21, TRA = "/")^2)) %>%  # Herfindahl Index
  fmutate(NHI = (HI - 1/N)/(1-1/N)) # Normalized Herfindahl Index


# Merging
ineq_96km2_hex_all <- IWI_ineq_96km2_hex %>% add_stub("IWI_", cols = -1L) %>% 
  merge(RWI_ineq_96km2_hex %>% add_stub("RWI_", cols = -1L), by = "cell", all = TRUE) %>% 
  merge(NL21_ineq_96km2_hex %>% add_stub("NL21_", cols = -1L), by = "cell", all = TRUE) %>% 
  fsubset(is.finite(IWI_mean) | is.finite(RWI_mean) | is.finite(NL21_mean))

fnobs(ineq_96km2_hex_all) # check if all matched, but seems so

# Checking correlations
ineq_96km2_hex_all %>% fselect(IWI_pop, RWI_pop) %>% na_omit() %$% cor(IWI_pop, RWI_pop)
ineq_96km2_hex_all %$% pwcor(IWI_w_GINI, RWI_w_GINI)

# Combined estimates 
set_collapse(na.rm = TRUE)
settransform(ineq_96km2_hex_all, 
   GINI_mean = psum(list(IWI_GINI, RWI_GINI, NL21_GINI) %>% TRA(proportions(1/fmean(.)), "*")),
   WGINI_mean = psum(list(IWI_w_GINI, RWI_w_GINI, NL21_w_GINI) %>% TRA(proportions(1/fmean(.)), "*")),
   TI_mean = psum(list(IWI_TI, RWI_TI, NL21_TI) %>% TRA(proportions(1/fmean(.)), "*")),
   HI_mean = psum(list(IWI_HI, RWI_HI, NL21_HI) %>% TRA(proportions(1/fmean(.)), "*")),
   NHI_mean = psum(list(IWI_NHI, RWI_NHI, NL21_NHI) %>% TRA(proportions(1/fmean(.)), "*"))
)

# Materialize grid cells
cells <- dgcellstogrid(world_96km2_hex, ineq_96km2_hex_all$cell) 
plot(cells) # plot them
identical(cells$seqnum, ineq_96km2_hex_all$cell) # Check that order is maintained

# Add them to the data file
ineq_96km2_hex_all %<>% add_vars(cells %>% fselect(-seqnum)) %>% 
  copyMostAttrib(cells) %>% st_make_valid()

# Saving as geo pakage
st_write(ineq_96km2_hex_all, "results/SA_ineq_96km2_hex_all.gpkg", delete_dsn = TRUE)

rm(cells, ineq_96km2_hex_all, IWI_ineq_96km2_hex, RWI_ineq_96km2_hex, NL21_ineq_96km2_hex, world_96km2_hex); gc()

#
# Inequality Within Grid at 1km Resolution Using Spatial Interpolation -----------------------------------------------------
#

fastverse_extend(dggridR, sf, s2, install = TRUE)
# High-res 1km grid
world_1km_hex <- dgconstruct(area = 1, metric = TRUE, resround = "nearest")
# I use Nightlights 2020 layer to remove cells outside the SA boundary (the layer was masked in GEE using a high-res country boundary, setting all values outside the mask to -9999)
NL20 <- terra::rast("data/south_africa_viirs_dnb_nightlights_v1_vcmslcfg_annual_median_composite/2020.tif") 
plot(NL20) # This is the grid we are looking at... 
NL20 %<>% as.data.frame(xy = TRUE) %>% fsubset(avg_rad %!=% -9999) %>% qDT() %>% frename(x = lon, y = lat)
# This gets the cell numbers
NL_1km_cells <- NL20 %$% dgGEO_to_SEQNUM(world_1km_hex, lon, lat)$seqnum
# This adds the cell centroids
NL_1km_cells <- dgSEQNUM_to_GEO(world_1km_hex, NL_1km_cells) %>% rm_stub("_deg", pre = FALSE) %>% 
                add_vars(cell = NL_1km_cells, pos = "front") %>% qDT()
rm(NL20); gc()

# Starting with IWI -----------------------------------------------------------------------------------------

# Takes around 2 min
system.time({
  IWI_1km_cells <- IWI_round %$% upsample_cells(world_1km_hex, lon, lat, 0.3, 50) %>% intersect(NL_1km_cells$cell)
})
IWI_1km_cells <- dgSEQNUM_to_GEO(world_1km_hex, IWI_1km_cells) %>% rm_stub("_deg", pre = FALSE) %>% 
                 add_vars(cell = IWI_1km_cells, pos = "front") %>% qDT()
  
# Now calculating IWI based on S2 geometries...
IWI_s2 = IWI_round %$% s2_lnglat(lon, lat)
IWI_1km_cells_s2 = IWI_1km_cells %$% s2_lnglat(lon, lat)

# This still takes quite long
s2_distance(IWI_s2, IWI_1km_cells_s2[[1]])

# Thus go degree-wise
# degree_grid = IWI_round %$% expand.grid(lon = floor(min(lon)):ceiling(max(lon)), 
#                                         lat = floor(min(lat)):ceiling(max(lat)))

degree_grid = IWI_round %$% expand.grid(lon = seq(floor(min(lon)), ceiling(max(lon)), 0.5),
                                        lat = seq(floor(min(lat)), ceiling(max(lat)), 0.5))
# Print out computation progresss
verbose = FALSE

IWI_1km_cells$GINI_5km = gridded_distance_interpol(degree_grid, IWI_round$IWI, IWI_s2, IWI_1km_cells_s2, 
                                               thresh = 5000, step = 0.5, tol = 0.1, FUN = w_gini, verbose = verbose)

IWI_1km_cells$WGINI_5km = gridded_distance_interpol(degree_grid, IWI_round$IWI, IWI_s2, IWI_1km_cells_s2, 
                                                thresh = 5000, step = 0.5, tol = 0.1, FUN = w_gini, w = IWI_round$pop, verbose = verbose)

IWI_1km_cells$GINI_10km = gridded_distance_interpol(degree_grid, IWI_round$IWI, IWI_s2, IWI_1km_cells_s2, 
                                               thresh = 10000, step = 0.5, tol = 0.1, FUN = w_gini, verbose = verbose)

IWI_1km_cells$WGINI_10km = gridded_distance_interpol(degree_grid, IWI_round$IWI, IWI_s2, IWI_1km_cells_s2, 
                                                thresh = 10000, step = 0.5, tol = 0.1, FUN = w_gini, w = IWI_round$pop, verbose = verbose)

IWI_1km_cells %<>% replace_Inf() %>% na_omit(cols = names(.) %like% "GINI", prop = 0.5)

IWI_1km_cells %>% gvr("GINI") %>% pwcor()

IWI_1km_cells %>% fwrite("results/SA_IWI_GINI_1km.csv")


# Same for RWI ---------------------------------------------------------------------------------------------------

system.time({
  RWI_1km_cells = RWI_round %$% upsample_cells(world_1km_hex, lon, lat, 0.3, 50) %>% intersect(NL_1km_cells$cell)
})
RWI_1km_cells <- dgSEQNUM_to_GEO(world_1km_hex, RWI_1km_cells) %>% rm_stub("_deg", pre = FALSE) %>% 
                 add_vars(cell = RWI_1km_cells, pos = "front") %>% qDT()

# Creating S2 geometries
RWI_s2 = RWI_round %$% s2_lnglat(lon, lat)
RWI_1km_cells_s2 = RWI_1km_cells %$% s2_lnglat(lon, lat)

# This still takes quite long
s2_distance(RWI_s2, RWI_1km_cells_s2[[1]])

degree_grid = RWI_round %$% expand.grid(lon = seq(floor(min(lon)), ceiling(max(lon)), 0.5),
                                        lat = seq(floor(min(lat)), ceiling(max(lat)), 0.5))

RWI_1km_cells$GINI_5km = gridded_distance_interpol(degree_grid, RWI_round$RWI, RWI_s2, RWI_1km_cells_s2, 
                                                thresh = 5000, step = 0.5, tol = 0.1, FUN = w_gini, verbose = verbose)

RWI_1km_cells$WGINI_5km = gridded_distance_interpol(degree_grid, RWI_round$RWI, RWI_s2, RWI_1km_cells_s2, 
                                                 thresh = 5000, step = 0.5, tol = 0.1, FUN = w_gini, w = RWI_round$pop, verbose = verbose)

RWI_1km_cells$GINI_10km = gridded_distance_interpol(degree_grid, RWI_round$RWI, RWI_s2, RWI_1km_cells_s2, 
                                                thresh = 10000, step = 0.5, tol = 0.1, FUN = w_gini, verbose = verbose)

RWI_1km_cells$WGINI_10km = gridded_distance_interpol(degree_grid, RWI_round$RWI, RWI_s2, RWI_1km_cells_s2, 
                                                 thresh = 10000, step = 0.5, tol = 0.1, FUN = w_gini, w = RWI_round$pop, verbose = verbose)

RWI_1km_cells %<>% replace_Inf() %>% na_omit(cols = names(.) %like% "GINI", prop = 0.5)

RWI_1km_cells %>% gvr("GINI") %>% pwcor()

RWI_1km_cells %>% fwrite("results/SA_RWI_GINI_1km.csv")

# Same for NL21 -------------------------------------------------------------------------------------------------------------------------

# Technically not necessary to upsample, but we want to take along zero nightlights cells between cells with positive nightlights
system.time({
  NL21_1km_cells = NL21_round %$% upsample_cells(world_1km_hex, lon, lat, 0.3, 10) %>% intersect(NL_1km_cells$cell)
})

NL21_1km_cells <- dgSEQNUM_to_GEO(world_1km_hex, NL21_1km_cells) %>% rm_stub("_deg", pre = FALSE) %>% 
                  add_vars(cell = NL21_1km_cells, pos = "front") %>% qDT()

# Creating S2 geometries
NL21_s2 = NL21_round %$% s2_lnglat(lon, lat)
NL21_1km_cells_s2 = NL21_1km_cells %$% s2_lnglat(lon, lat)

# This takes quite long
s2_distance(NL21_s2, NL21_1km_cells_s2[[1]])

degree_grid = NL21_round %$% expand.grid(lon = seq(floor(min(lon)), ceiling(max(lon)), 0.25),
                                         lat = seq(floor(min(lat)), ceiling(max(lat)), 0.25))

# These estimates can take a few hours to complete, even on a fast machine... 
NL21_1km_cells$GINI_5km = gridded_distance_interpol(degree_grid, NL21_round$NL21, NL21_s2, NL21_1km_cells_s2, 
                                               thresh = 5000, step = 0.25, tol = 0.05, FUN = w_gini, verbose = TRUE)

NL21_1km_cells$WGINI_5km = gridded_distance_interpol(degree_grid, NL21_round$NL21, NL21_s2, NL21_1km_cells_s2, 
                                               thresh = 5000, step = 0.25, tol = 0.05, FUN = w_gini, w = NL21_round$pop, verbose = TRUE)

NL21_1km_cells$GINI_10km = gridded_distance_interpol(degree_grid, NL21_round$NL21, NL21_s2, NL21_1km_cells_s2, 
                                               thresh = 10000, step = 0.25, tol = 0.1, FUN = w_gini, verbose = TRUE)

NL21_1km_cells$WGINI_10km = gridded_distance_interpol(degree_grid, NL21_round$NL21, NL21_s2, NL21_1km_cells_s2, 
                                               thresh = 10000, step = 0.25, tol = 0.1, FUN = w_gini, w = NL21_round$pop, verbose = TRUE)

NL21_1km_cells %<>% replace_Inf() %>% na_omit(cols = names(.) %like% "GINI", prop = 0.5)

NL21_1km_cells %>% gvr("GINI") %>% pwcor()

NL21_1km_cells %>% fwrite("results/SA_NL21_GINI_1km.csv")
