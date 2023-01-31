###########################################
# Population Weighted Geospatial Inequality
###########################################

library(fastverse)
fastverse_extend(qs) 
set_collapse(nthreads = 4, na.rm = FALSE)

gini_wiki <- function(x) {1 + 2/(length(x)-1) * (sum(seq_along(x)*sort(x)) / sum(x) - length(x))} 
gini_noss <- function(x) 2/length(x) * sum(seq_along(x)*sort(x)) / sum(x) - (length(x)+1)/length(x) # same as in Galimberti et al. (2020)

# S <- function(p, q) {q/2 * (2*p + q - 1)}
# all.equal(sum(30:100), S(30, 70+1))
# Sn <- function(n) n*(n+1)/2
# all.equal(sum(1:100), Sn(100))
# Skn <- function(k, n) (n*(n+1) - (k-1)*k)/2 
# all.equal(Skn(30, 100), sum(30:100))
# Skpq <- function(k, q) q/2 * (2*k + q + 1) + k
# all.equal(Skpq(30, 70), sum(30:100))
Skp1qm1 <- function(k, q) (q-1)/2 * (2*(k+1) + q) + k + 1
all.equal(Skp1qm1(30-1, 70+1), sum(30:100))

w_gini <- function(x, w, sscor = FALSE) {
  si = .Internal(qsort(x, TRUE))
  w = w[si$ix]
  x = si$x
  sw = sum(w)
  csw = cumsum(w)
  sx = Skp1qm1(c(0, csw[-length(csw)]), w) # Skpq(c(0, csw[-length(csw)])+1, w-1)
  if(sscor) return(1 + 2/(sw-1)*(sum(sx*x) / sum(x*w) - sw)) # TODO: see paper again, correction should be different 
  2/sw * sum(sx*x) / sum(x*w) - (sw+1)/sw
}


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
