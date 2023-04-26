
library(fastverse)
fastverse_extend(sf, install = TRUE)


# 96km2 Hexagonal Grid -------------------------------------------------------------------------

ineq_96km2_hex_all <- st_read("results/SA_ineq_96km2_hex_all.gpkg")

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


# 1km Interpolations ---------------------------------------------------------------------

GINI_1km <- fread("results/SA_IWI_GINI_1km.csv") %>% add_stub("IWI_", cols = -(1:3)) %>%
  merge(fread("results/SA_RWI_GINI_1km.csv") %>% slt(-lon, -lat) %>% add_stub("RWI_", cols = -1L), by = "cell") %>%
  merge(fread("results/SA_NL21_GINI_1km.csv") %>% slt(-lon, -lat) %>% add_stub("NL_", cols = -1L), by = "cell") %>%
  ftransform(GINI_5km_mean = psum(list(IWI_GINI_5km, RWI_GINI_5km, NL_GINI_5km) %>% TRA(proportions(1/fmean(.)), "*")),
             WGINI_5km_mean = psum(list(IWI_WGINI_5km, RWI_WGINI_5km, NL_WGINI_5km) %>% TRA(proportions(1/fmean(.)), "*")), 
             GINI_10km_mean = psum(list(IWI_GINI_5km, RWI_GINI_5km, NL_GINI_5km) %>% TRA(proportions(1/fmean(.)), "*")),
             WGINI_10km_mean = psum(list(IWI_WGINI_5km, RWI_WGINI_5km, NL_WGINI_5km) %>% TRA(proportions(1/fmean(.)), "*")))

GINI_1km %>% gvr("_GINI|^GINI_") %>% pwcor()
GINI_1km %>% gvr("_WGINI|^WGINI_") %>% pwcor()
GINI_1km %>% gvr("_mean") %>% pwcor()

# Only 5km
GINI_1km %>% gvr("_WGINI_5km|^WGINI_5km") %>% pwcor(N = TRUE, P = TRUE)


# Relating 96km2 HEX GINI's to STP3 ---------------------------------------------------------

H7_grid_all <- st_read("results/STP3_Hex7_grid.gpkg")

# Cape Town 96km2 hexagons
ind_ct <- st_contains(ineq_96km2_hex_all$geom, st_centroid(subset(H7_grid_all, CAT_B == "CPT", geom))) %>% 
  lengths() %>% is_greater_than(0) %>% which()

# Matching with Uber Hexagons
ind <- st_contains(ineq_96km2_hex_all$geom, st_centroid(H7_grid_all$geom)) 
H7_grid_all$HEX96_row <- NA_integer_
for(i in seq_along(ind)) setv(H7_grid_all$HEX96_row, ind[[i]], i, vind1 = TRUE)
descr(H7_grid_all$HEX96_row) # 1% missing

# Adding STP3 data to 96km2 hexagons
tmp <- H7_grid_all %>% qDT() %>% 
  fsubset(!is.na(HEX96_row), MedianIncome_2020, FTE_2020, gini_2020, HEX96_row) %>% 
  fgroup_by(HEX96_row) %>% fmean(FTE_2020)

ineq_96km2_hex_all[tmp$HEX96_row, names(tmp)[-1]] <- tmp[, -1]
rm(tmp)

# Plotting various weighted GINI's in Cape Town (again wide graphics device to have them all in one row)
plot(ineq_96km2_hex_all[ind_ct, .c(gini_2020, IWI_w_GINI, RWI_w_GINI, NL21_w_GINI)] %>% 
     frename(gini_2020 = STP3_w_GINI_2020), lwd = 0.1)

dev.copy(pdf, "figures/Hex95km_CompGini.pdf", width = 11.69, height = 5)
dev.off()

# Correlations
ineq_96km2_hex_all %>% qDT() %>% gvr("w_GINI|gini_2020|WGINI") %>%
  frename(gini_2020 = STP3_w_GINI_2020) %>%
  pwcor(N = TRUE, P = TRUE) %>% 
  print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE) %>%
  xtable::xtable() %>% print(booktabs = TRUE)


# Relating 1km2 Interpolated GINI's to STP3 ---------------------------------------------------------

GINI_1km_sf <- GINI_1km %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Matching to Uber Hexagons
ind <- st_contains(H7_grid_all$geom, GINI_1km_sf$geometry)
GINI_1km_sf$UberH7 <- NA_integer_
for(i in seq_along(ind)) setv(GINI_1km_sf$UberH7, ind[[i]], i, vind1 = TRUE)
descr(GINI_1km_sf$UberH7) # No missing values !
names(GINI_1km_sf)

# Adding interpolated GINI's to Uber hexagons
tmp <- GINI_1km_sf %>% qDT() %>% 
  fsubset(!is.na(UberH7), -cell) %>% 
  fgroup_by(UberH7) %>% num_vars() %>% fmean()

H7_grid_all[tmp$UberH7, names(tmp)[-1]] <- tmp[, -1]
names(H7_grid_all)

# Plotting various weighted GINI's in Cape Town (again wide graphics device to have them all in one row)
H7_grid_all[H7_grid_all$CAT_B == "CPT", ] %>% 
  gvr("w_GINI|gini_2020|WGINI_5km") %>% 
  fselect(-WGINI_5km_mean) %>% plot(lwd = 1)

dev.copy(pdf, "figures/Hex1_3km_CompGini.pdf", width = 11.69, height = 5)
dev.off()

# Correlations
H7_grid_all %>% qDT() %>% gvr("w_GINI|gini_2020|WGINI_5km") %>%
  colorder(gini_2020, pos = "end") %>% 
  frename(gini_2020 = STP3_w_GINI_2020) %>%
  pwcor(N = TRUE, P = TRUE) %>% 
  print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE) %>%
  xtable::xtable() %>% print(booktabs = TRUE)
