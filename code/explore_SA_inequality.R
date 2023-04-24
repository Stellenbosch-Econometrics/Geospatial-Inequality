
library(fastverse)
fastverse_extend(sf, install = TRUE)

#
# Comparison of Estimates ---------------------------------------------------------------
# 

# 96km2 hex grid -------------------------------------------------------------------------
ineq_96km2_hex_all <- st_read("results/hex_grids/SA_ineq_96km2_hex_all.gpkg")

ineq_96km2_hex_all %>% qDT() %>% gvr("pop") %>% pwcor()

ineq_96km2_hex_all %>% qDT() %>% gvr("mean") %>% gvr("w", invert = TRUE) %>% pwcor()
ineq_96km2_hex_all %>% qDT() %>% gvr("w_mean") %>% pwcor()

ineq_96km2_hex_all %>% qDT() %>% gvr("median") %>% gvr("w", invert = TRUE) %>% pwcor()
ineq_96km2_hex_all %>% qDT() %>% gvr("w_median") %>% pwcor()

ineq_96km2_hex_all %>% qDT() %>% gvr("GINI") %>% gvr("w|WGINI", invert = TRUE) %>% pwcor()
ineq_96km2_hex_all %>% qDT() %>% gvr("w_GINI|WGINI") %>% pwcor()

ineq_96km2_hex_all %>% qDT() %>% gvr("TI") %>% pwcor()
ineq_96km2_hex_all %>% qDT() %>% gvr("NHI") %>% pwcor()

ineq_96km2_hex_all %>% qDT() %>% fselect(GINI_mean, WGINI_mean, TI_mean, HI_mean, NHI_mean) %>% pwcor()


# 1km interpolations ---------------------------------------------------------------------

GINI_1km <- fread("results/1km_interpolations/IWI_GINI_1km_hex.csv") %>% add_stub("IWI_", cols = 4:5) %>%
  merge(fread("results/1km_interpolations/RWI_GINI_1km_hex.csv") %>% slt(-lon_deg, -lat_deg) %>% add_stub("RWI_", cols = -1L), by = "cell") %>%
  merge(fread("results/1km_interpolations/NL20_GINI_1km_hex.csv") %>% slt(-lon_deg, -lat_deg) %>% add_stub("NL_", cols = -1L), by = "cell") %>%
  ftransform(GINI_mean = psum(list(IWI_GINI, RWI_GINI, NL_GINI) %>% TRA(proportions(1/fmean(.)), "*")),
             WGINI_mean = psum(list(IWI_WGINI, RWI_WGINI, NL_WGINI) %>% TRA(proportions(1/fmean(.)), "*")))

GINI_1km %>% gvr("_GINI|^GINI_") %>% pwcor()
GINI_1km %>% gvr("_WGINI|^WGINI_") %>% pwcor()
GINI_1kms %>% gvr("_mean") %>% pwcor()

#
### Exploration -----------------------------------------
#


#################################################################
# This script relates interpolation based GINI's
# to the STP3, to produce some tables and figures in the slides 
#################################################################

# Cape town 
ind_ct = st_contains(SA_ineq_96km2_hex_all$geom, st_centroid(subset(H7_grid_all, CAT_B == "CPT", geom))) %>% 
  lengths() %>% is_greater_than(0) %>% which()

# Matching with Uber Hexagons
ind = st_contains(SA_ineq_96km2_hex_all$geom, st_centroid(H7_grid_all$geom)) 
H7_grid_all$HEX96_row <- NA_integer_
for(i in seq_along(ind)) setv(H7_grid_all$HEX96_row, ind[[i]], i, vind1 = TRUE)

names(H7_grid_all)
set_collapse(sort = TRUE, na.rm = TRUE)

tmp =  H7_grid_all %>% qDT() %>% 
  fsubset(!is.na(HEX96_row), MedianIncome_2020, FTE_2020, gini_2020, HEX96_row) %>% 
  fgroup_by(HEX96_row) %>% fmean(FTE_2020)

SA_ineq_96km2_hex_all[tmp$HEX96_row, names(tmp)[-1]] <- tmp[, -1]


ind_all = ind %>% lengths() %>% is_greater_than(0) %>% which()

SA_ineq_96km2_hex_all <- SA_ineq_96km2_hex_all[ind_all, ]

plot(SA_ineq_96km2_hex_all[ind_ct, .c(gini_2020, IWI_w_GINI, RWI_w_GINI, NL21_w_GINI)] %>% 
     frename(gini_2020 = STP3_w_GINI_2020), lwd = 0.1)

dev.copy(pdf, "figures/Hex95km_CompGini.pdf", width = 11.69, height = 5)
dev.off()


SA_ineq_96km2_hex_all %>% qDT() %>% gvr("w_GINI|gini_2020|WGINI") %>%
  frename(gini_2020 = STP3_w_GINI_2020) %>%
  pwcor(N = TRUE, P = TRUE) %>% 
  print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE) %>%
  xtable::xtable() %>% print(booktabs = TRUE)


# Interpolations -------------------------------------

SA_1km_5km_radius <- fread("data/1km_interpolations/SA_GINI_1km_hex_5km_radius.csv")

SA_1km_5km_radius %>% gvr("_GINI|^GINI_") %>% pwcor()
SA_1km_5km_radius %>% gvr("_WGINI|^WGINI_") %>% pwcor()
SA_1km_5km_radius %>% gvr("_mean") %>% pwcor()

SA_1km_5km_radius %>% gvr("_WGINI|^WGINI_") %>% 
  pwcor(N = TRUE, P = TRUE) %>% 
  print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE) %>%
  xtable::xtable() %>% print(booktabs = TRUE)


# Matching with Uber Hexagons -----------------------

SA_1km_5km_radius_sf <- SA_1km_5km_radius %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

ind = st_contains(H7_grid_all$geom, SA_1km_5km_radius_sf$geometry)

SA_1km_5km_radius_sf$Uberh7 <- NA_integer_
for(i in seq_along(ind)) setv(SA_1km_5km_radius_sf$Uberh7, ind[[i]], i, vind1 = TRUE)

names(SA_1km_5km_radius_sf)
set_collapse(sort = TRUE, na.rm = TRUE)

tmp = SA_1km_5km_radius_sf %>% qDT() %>% 
  fsubset(!is.na(Uberh7)) %>% 
  fgroup_by(Uberh7) %>% num_vars() %>% fmean()

H7_grid_all[tmp$Uberh7, names(tmp)[-1]] <- tmp[, -1]

names(H7_grid_all)

H7_grid_all[H7_grid_all$CAT_B == "CPT", ] %>% gvr("w_GINI|gini_2020|WGINI") %>% 
  fselect(-WGINI_mean)  %>% plot(lwd = 1)

dev.copy(pdf, "figures/Hex1_3km_CompGini.pdf", width = 11.69, height = 5)
dev.off()


H7_grid_all %>% qDT() %>% gvr("w_GINI|gini_2020|WGINI") %>%
  colorder(gini_2020, pos = "end") %>% 
  frename(gini_2020 = STP3_w_GINI_2020) %>%
  pwcor(N = TRUE, P = TRUE) %>% 
  print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE) %>%
  xtable::xtable() %>% print(booktabs = TRUE)
