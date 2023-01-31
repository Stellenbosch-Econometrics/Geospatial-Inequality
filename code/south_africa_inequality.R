############################
# Analysis for South Africa
############################

library(fastverse)
set_collapse(nthreads = 4, na.rm = FALSE)

SAIWI <- fread("data/South Africa_estimated_wealth_index.csv") %>% frename(estimated_IWI = IWI)
SARWI <- fread("data/South Africa_RWI.csv") %>% frename(longitude = lon, latitude = lat, rwi = RWI)

source("code/gini_helper.R")
source("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/code/osm_helpers.R")

SARWI %$% round_to_kms_exact(lon, lat, 2.4) %>% GRPN(FALSE) %>% descr()
SAIWI %$% round_to_kms_exact(lon, lat, 1.7) %>% GRPN(FALSE) %>% descr()

SARWI_round <- SARWI %>% ftransform(round_to_kms_exact(lon, lat, 2.4)) %>% 
  fgroup_by(lon, lat, sort = FALSE) %>% fselect(RWI) %>% fmean() %>% 
  merge(WPOP %>% ftransform(round_to_kms_exact(lon, lat, 2.4)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat))

SAIWI_round <- SAIWI %>% ftransform(round_to_kms_exact(lon, lat, 1.7)) %>% 
  fgroup_by(lon, lat, sort = FALSE) %>% fselect(IWI) %>% fmean() %>% 
  merge(WPOP %>% ftransform(round_to_kms_exact(lon, lat, 1.7)) %>% 
          fgroup_by(lon, lat, sort = FALSE) %>% fsum(), by = .c(lon, lat))

settransform(SARWI_round, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)
settransform(SAIWI_round, cell = dgGEO_to_SEQNUM(world_10km_hex, lon, lat)$seqnum)

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


# Merging
SA_ineq_10km_hex_all <- SAIWI_ineq_10km_hex %>% add_stub("IWI_", cols = -1L) %>% 
  merge(SARWI_ineq_10km_hex %>% add_stub("RWI_", cols = -1L), by = "cell", all = TRUE) %>% 
  fsubset(is.finite(IWI_mean) | is.finite(RWI_mean))

fnobs(SA_ineq_10km_hex_all) # check if all matched, but seems so
SA_ineq_10km_hex_all %>% fselect(IWI_pop, RWI_pop) %>% na_omit() %$% cor(IWI_pop, RWI_pop)

# Checking correlations
SA_ineq_10km_hex_all %$% pwcor(IWI_w_GINI, RWI_w_GINI)

# Adding coordinates
add_vars(SA_ineq_10km_hex_all) <- fsubset(africa_10km_hex_sf, ckmatch(SA_ineq_10km_hex_all$cell, cell), -cell)
SA_ineq_10km_hex_all %<>% st_as_sf()

# Saving as geo pakage
st_write(SA_ineq_10km_hex_all, "data/SA_ineq_10km_hex_all.gpkg")