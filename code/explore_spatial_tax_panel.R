############################################################
# Explore Inequality Estimates on the Spatial Tax Panel v3.7
############################################################

library(fastverse)
fastverse_extend(sf, install = TRUE)

range_title <- function(title, x, digits = 2) {
  paste(title, do.call(sprintf, as.list(c("[%s, %s]", round(.range(x), digits)))))
}

#
### Municipal Level ----------------------------------------------------------------------
#

MunGeo_data <- st_read("results/STP3_MunGeo_data.gpkg")

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

# Plotting wealth estimates 
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

# Plotting Raw GINI estimates
oldpar <- par(mfrow = c(2, 3), oma = c(0,0,0,0), mar = c(0,0,0,0))
plot(MunGeo_data["IWI_GINI"], main = range_title("International Wealth Index GINI", MunGeo_data$IWI_GINI), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["RWI_GINI"], main = range_title("Relative Wealth Index GINI", MunGeo_data$RWI_GINI), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["NL20_GINI"], main = range_title("VIIRS-DNB 2020 GINI", MunGeo_data$NL20_GINI), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["IWI_W_GINI"], main = range_title("International Wealth Index Weighted GINI", MunGeo_data$IWI_W_GINI), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["RWI_W_GINI"], main = range_title("Relative Wealth Index Weighted GINI", MunGeo_data$RWI_W_GINI), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["NL20_W_GINI"], main = range_title("VIIRS-DNB 2020 Weighted GINI", MunGeo_data$NL20_W_GINI), border = NA, key.pos = NULL, reset = FALSE)
par(oldpar)

dev.copy(pdf, "figures/STP3_municipal_gini_ests_raw.pdf", width = 11.69, height = 8.27)
dev.off()

# Only examining the most populated municipalities
populated_mun <- MunGeo_data %$% CAT_B[POP20 > 200000]

# Correlations of Raw GINI estimates
MunGeo_data %>% qDT() %>% 
  gvr("CAT_B|GINI$|gini_avg") %>% 
  fsubset(CAT_B %in% populated_mun) %>% # comment out this line for all municipalities
  num_vars() %>%
  colorderv(neworder = "W_", regex = TRUE) %>%
  colorder(gini_avg, pos = "end") %>% 
  pwcor(P = TRUE) %>% 
  print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE) %>%
  xtable::xtable() %>% print(booktabs = TRUE)
  
# Plotting Calibrated GINI estimate
oldpar <- par(mfrow = c(2, 3), oma = c(0,0,0,0), mar = c(0,0,0,0))
plot(MunGeo_data["IWI_GINI_calib"], main = range_title("International Wealth Index GINI", MunGeo_data$IWI_GINI_calib), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["RWI_GINI_calib"], main = range_title("Relative Wealth Index GINI", MunGeo_data$RWI_GINI_calib), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["NL20_GINI_calib"], main = range_title("VIIRS-DNB 2020 GINI", MunGeo_data$NL20_GINI_calib), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["IWI_W_GINI_calib"], main = range_title("International Wealth Index Weighted GINI", MunGeo_data$IWI_W_GINI_calib), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["RWI_W_GINI_calib"], main = range_title("Relative Wealth Index Weighted GINI", MunGeo_data$RWI_W_GINI_calib), border = NA, key.pos = NULL, reset = FALSE)
plot(MunGeo_data["NL20_W_GINI_calib"], main = range_title("VIIRS-DNB 2020 Weighted GINI", MunGeo_data$NL20_W_GINI_calib), border = NA, key.pos = NULL, reset = FALSE)
par(oldpar)

dev.copy(pdf, "figures/STP3_municipal_gini_ests_calib.pdf", width = 11.69, height = 8.27)
dev.off()

# Correlations of Calibrated GINI estimates
MunGeo_data %>% qDT() %>% 
  gvr("CAT_B|GINI_calib$|gini_avg") %>% 
  fsubset(CAT_B %in% populated_mun) %>% # comment out this line for all municipalities
  num_vars() %>% 
  colorderv(neworder = "W_", regex = TRUE) %>% 
  colorder(gini_avg, pos = "end") %>% 
  pwcor(P = TRUE) %>% 
  print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE) %>%
  xtable::xtable() %>% print(booktabs = TRUE)



#
### 5.5km2 Uber Hexagon Level --------------------------------------------------------------
#

H7_grid_all <- st_read("results/STP3_Hex7_grid.gpkg")

# Explore correlations
H7_grid_all %>% qDT() %>% 
  fselect(IWI, RWI, NL20, POP20, MedianIncome_2020, FTE_2020, gini_2020) %>% 
  pwcor(P = TRUE, N = TRUE) %>% print(digits = 3, sig.level = 0.05, show = "lower.tri")

H7_grid_all %>% qDT() %>% 
  fmutate(lat_lon_int = finteraction(as.integer(lat), as.integer(lon))) %>% 
  fcomputev(c(IWI, RWI, NL20, MedianIncome_2020, FTE_2020, gini_2020), 
            fscale, lat_lon_int, na.rm = TRUE) %>% 
  pwcor(P = TRUE, N = TRUE) %>% print(digits = 3, sig.level = 0.05, show = "lower.tri")

H7_grid_all %>% qDT() %>% 
  fcomputev(c(IWI, RWI, NL20, MedianIncome_2020, FTE_2020, gini_2020), 
            fscale, CAT_B, na.rm = TRUE) %>% 
  pwcor(P = TRUE, N = TRUE) %>% print(digits = 3, sig.level = 0.05, show = "lower.tri")

# Plotting 
H7_grid_all %>% subset(CAT_B == "CPT", MedianIncome_2020) %>% plot()

# The graphics device needs to be wide and narrow for all graphs to go in one row
H7_grid_all %>% subset(CAT_B == "CPT", c(MedianIncome_2020, IWI, RWI, NL20)) %>% plot(lwd = 0.1)

dev.copy(pdf, "figures/STP3_CompInc_UberHex.pdf", width = 11.69, height = 5)
dev.off()

# Export Correlations: overall and within municipalities
H7_grid_all %>% qDT() %>% 
  fselect(MedInc20 = MedianIncome_2020, IWI, RWI, NL20, CAT_B) %>% {
    cbind(Overall = num_vars(.) %>% pwcor(N = TRUE, P = TRUE) %>% 
            print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE),
          Within = STD(., ~CAT_B, stub = FALSE, keep.by = FALSE) %>% pwcor(N = TRUE, P = TRUE) %>% 
            print(digits = 3, sig.level = 0.05, show = "lower.tri", return = TRUE))
  } %>% xtable::xtable() %>% print(booktabs = TRUE)


