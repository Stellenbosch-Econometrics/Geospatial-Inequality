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

#
# 10km Hexagonal grid -----------------------------------------------
#
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

#
# Now: Grid at 1km resolution ----------------------------------------------------------------
#
fastverse_extend(qs, dggridR, sf, s2)
world_1km_hex <- dgconstruct(area = 1, metric = TRUE, resround = "nearest")
# africa_10km_hex_sf %>% subset(cell %in% funique(c(SAIWI_round$cell, SARWI_round$cell))) %>% 
#   st_write("data/temporary/south_africa_IWI_RWI_hex_10km_cells.shp")
# centroids = africa_10km_hex_centroids %>% subset(cell %in% funique(c(SAIWI_round$cell, SARWI_round$cell))) %>% 
#   st_coordinates() %>% qDF()
# SA_1km_hex_sp <- dgshptogrid(world_1km_hex, "data/temporary/south_africa_IWI_hex_10km_cells.shp", cellsize = 0.01, frame = FALSE)
# SA_1km_hex_sp <- st_as_sf(africa_10km_hex_sp)
# rm(SA_1km_hex_sp); gc()
# # Convert to WGS84
# africa_10km_hex_sf %<>% st_transform(4326)
# # Trying to upsample with sf: Too expensive
# upsample_centroids = 
#   st_make_grid(x = africa_10km_hex_sf %>% subset(cell %in% funique(c(SAIWI_round$cell, SARWI_round$cell))),
#                cellsize = 0.005, what = "centers", square = FALSE)


# Thus need to do things manually...
upsample_cells <- function(dggrid, lon, lat, km, times) {
  degrees = km / (40075.017 / 360)  # Gets the degree-distance of the kms at the equator
  cells = dgGEO_to_SEQNUM(dggrid, lon, lat)$seqnum
  init = -degrees * times / 2
  lat2 = lat + init; lon2 = lon + init
  for (j in seq_len(times)) {
    for (k in seq_len(times)) {
      lon2 %+=% degrees
      cells = funique(c(cells, dgGEO_to_SEQNUM(dggrid, lon2, lat2)$seqnum))
    }
    lat2 %+=% degrees
    lon2 = lon + init
    cells = funique(c(cells, dgGEO_to_SEQNUM(dggrid, lon2, lat2)$seqnum))
  }
  return(cells)
}

# Takes around 1 min each...
system.time({
  SAIWI_cells = SAIWI_round %$% upsample_cells(world_1km_hex, lon, lat, 0.3, 50)
})
SAIWI_cells = add_vars(qDT(dgSEQNUM_to_GEO(world_1km_hex, SAIWI_cells)), cell = SAIWI_cells, pos = "front") %>% rm_stub("_deg", FALSE)

system.time({
  SARWI_cells = SARWI_round %$% upsample_cells(world_1km_hex, lon, lat, 0.3, 50)
})
SARWI_cells = add_vars(qDT(dgSEQNUM_to_GEO(world_1km_hex, SARWI_cells)), cell = SARWI_cells, pos = "front")

# Now calculating based on S2 geometries...

SAIWI_s2 = SAIWI_round %$% s2_lnglat(lon, lat)
SAIWI_cells_s2 = SAIWI_cells %$% s2_lnglat(lon_deg, lat_deg)

# This still takes quite long
s2_distance(SAIWI_s2, SAIWI_cells_s2[[1]])
# Thus go degree-wise
degree_grid = SAIWI_round %$% expand.grid(lon = floor(min(lon)):ceiling(max(lon)), 
                                          lat = floor(min(lat)):ceiling(max(lat)))

gridded_distance_calc <- function(degree_grid, y, points_s2, cells_s2, FUN, thresh, step = 1L, tol = 0.1, w = NULL, ..., verbose = FALSE) {
  if(length(y) != length(points_s2)) stop("length(y) must match length(points_s2)")
  degrees = qM(degree_grid)
  if(!identical(colnames(degrees), c("lon", "lat"))) stop("degree_grid must be a data frame with columns named 'lon' and 'lat'")
  wnull = is.null(w)
  if(!wnull && length(w) != length(y)) stop("length(w) must match length(y)")

  res = numeric(length(y))
  px = s2_x(points_s2); py = s2_y(points_s2); cx = s2_x(cells_s2); cy = s2_y(cells_s2)
  
  for(d in mrtl(degrees)) {
      if(verbose) cat("lon =", d[1L], "lat =", d[2L], fill = TRUE)
      which_cells = which(between(cx, d[1L], d[1L]+step) & between(cy, d[2L], d[2L]+step))
      which_points = which(between(px, d[1L]-tol, d[1L]+(tol+step)) & between(py, d[2L]-tol, d[2L]+(tol+step)))
      points_s2_subset = points_s2[which_points]
      y_subset = y[which_points]
      if(wnull) {
        for(i in which_cells) {
          disti = s2_distance(points_s2_subset, cells_s2[i])
          res[i] = FUN(y_subset[disti < thresh], ...)
        }
      } else {
        w_subset = w[which_points]
        for(i in which_cells) {
          disti = which(s2_distance(points_s2_subset, cells_s2[i]) < thresh)
          res[i] = FUN(y_subset[disti], w_subset[disti], ...)
        }
      }
  }
  return(res)
}

degree_grid = SAIWI_round %$% expand.grid(lon = seq(floor(min(lon)), ceiling(max(lon)), 0.5), 
                                          lat = seq(floor(min(lat)), ceiling(max(lat)), 0.5))

SAIWI_cells$GINI = gridded_distance_calc(degree_grid, SAIWI_round$IWI, SAIWI_s2, SAIWI_cells_s2, 
                                   thresh = 10000, step = 0.5, tol = 0.1, FUN = gini_wiki, verbose = TRUE)

SAIWI_cells$WGINI = gridded_distance_calc(degree_grid, SAIWI_round$IWI, SAIWI_s2, SAIWI_cells_s2, 
                                   thresh = 10000, step = 0.5, tol = 0.1, FUN = w_gini, w = SAIWI_round$pop, verbose = TRUE)

fwrite(SAIWI_cells, "data/SA_IWI_GINI_1km_hex.csv")

