###########################################
# South Africa Nightlights GINI Time Series
###########################################

library(fastverse)
set_collapse(nthreads = 4, na.rm = FALSE, sort = FALSE)
source("code/gini.R")
source("code/spatial_helpers.R")

nl_annual_path <- "data/south_africa_viirs_dnb_nightlights_v1_vcmslcfg_annual_median_composite"
nl_files <- list.files(nl_annual_path)

nl_data <- nl_files %>% set_names(substr(., 1, 4)) %>% 
  lapply(function(x) paste0(nl_annual_path, "/", x) %>% 
           terra::rast() %>% 
           terra::as.data.frame(xy = TRUE) %>% 
           fsubset(avg_rad %!=% -9999)) %>% 
  unlist2d("year", id.factor = TRUE, DT = TRUE) %>%
  frename(x = lon, y = lat) %>% 
  fmutate(year = as.integer(levels(year))[year])

# What population to measure? Daytime or nighttime?
pop_path <- "data/WPOP_SA_1km_UNadj"
pop_files <- list.files(pop_path)

pop_data <- pop_files %>% set_names(substr(., 9, 12)) %>% 
  lapply(function(x) paste0(pop_path, "/", x) %>% 
           terra::rast() %>% 
           terra::as.data.frame(xy = TRUE) %>% 
           set_names(c("lon", "lat", "pop"))) %>% 
  unlist2d("year", id.factor = TRUE, DT = TRUE) %>% 
  fmutate(year = as.integer(levels(year))[year])

# Forecasting Population 

pop_forecast <- pop_data %>% 
  roworder(year) %>% 
  fgroup_by(lat, lon) %>% 
  fmutate(dm_year = fwithin(year)) %>% 
  fsummarise(pop_2020 = flast(pop),
             beta = fsum(pop, dm_year) %/=% fsum(dm_year^2)) %>% 
  fmutate(pop_2021 = pop_2020 + beta, 
          pop_2022 = pop_2021 + beta, 
          beta = NULL)

head(pop_forecast)
# COVID cause migration into big cities and out of rural cities...

pop_data_forecast <- rbind(pop_data,
  pop_forecast %>% fselect(-pop_2020) %>% rm_stub("pop_") %>% 
  melt(1:2, variable.name = "year", value.name = "pop") %>% 
  fmutate(year = as.integer(levels(year))[year]) %>% 
  colorder(year))

# Checks for Spatial Matching

pop_data %>% fselect(year, lat, lon) %>% any_duplicated()
nl_data %>% fselect(year, lat, lon) %>% any_duplicated()

pop_data %>% ftransform(round_to_kms_fast(lon, lat, 0.63)) %>% 
  fselect(year, lat, lon) %>% any_duplicated()

nl_data %>% ftransform(round_to_kms_fast(lon, lat, 0.63)) %>% 
  fselect(year, lat, lon) %>% fduplicated() %>% qtable()

# Spatial Matching
nl_pop_data <- pop_data_forecast %>% 
  ftransform(round_to_kms_fast(lon, lat, 0.63)) %>% 
  merge(nl_data %>% ftransform(round_to_kms_fast(lon, lat, 0.63)) %>% 
         fgroup_by(year, lat, lon) %>% fmean(), by = .c(year, lat, lon)) %>% 
  roworder(year, lat, lon)

qsu(nl_pop_data)
rm(pop_data, pop_forecast, pop_data_forecast, nl_data)
gc()

qs::qsave(nl_pop_data, "data/sa_nl_pop_ts_data.qs")

#
### Now GINI Calculations -------------------------------
#

nl_pop_data <- qs::qread("data/sa_nl_pop_ts_data.qs")

set_collapse(sort = TRUE)

raw_gini <- nl_pop_data %>% 
  fsubset(pop > 0 & avg_rad > 0) %>%
  # ftransformv(c(pop, avg_rad), fmin, TRA = "-") %>%
  fgroup_by(year) %>% 
  fsummarise(gini = gini_noss(avg_rad)*100, 
             w_gini = w_gini(avg_rad, pop)*100) 

# raw_gini %$% ts(cbind(gini, w_gini), start = year[1]) %>% plot()

library(ggplot2)
raw_gini %>% melt("year") %>% 
  frename(year = Year, value = Value, variable = Series) %>% 
  ggplot(aes(x = Year, y = Value, colour = Series)) + geom_line() +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  scale_color_brewer(palette = "Paired") +
  theme_bw(base_size = 14) + 
  ggtitle("Uncalibrated VIIRS-DNB Nightlights Based GINI Estimates")

dev.copy(pdf, "figures/Raw_VIIRS_SA_GINI_TimeSeries.pdf", width = 9.27, height = 5.83)
dev.off()

#
### Inequality Estimates --------------------------------
#

# World Inequality Database
SA_WID <- fread("data/SA_GINI/WID/WID_SA_GINI.csv") %T>% with(value %*=% 100)

# Income based: 0.75 since 2014. Available on the website but not in the downloaded data...
library(ggplot2) # Problem: Wealth based GINI above 1? -> Possible with negative wealth....
SA_WID %>% ggplot(aes(x = year, y = value, colour = variable)) + geom_line() + 
  scale_y_continuous(limits = c(50, 110))

# World Bank Estimate
SA_WB <- africamonitor::am_data(series = "SI_POV_GINI", ctry = "ZAF") 
SA_WB %>% ggplot(aes(x = Date, y = SI_POV_GINI)) + geom_line() + scale_y_continuous(limits = c(50, 110))

# Standardized World Inequality Database: Gini based on disposable and Market Income
SA_SWID <- fread("data/SA_GINI/SWID/swiid9_4/swiid9_4_summary.csv") %>% 
           fsubset(country == "South Africa", year, gini_disp, gini_mkt) 

SA_SWID %>% melt("year") %>% ggplot(aes(x = year, y = value, colour = variable)) + geom_line() + 
  scale_y_continuous(limits = c(50, 110))
# -> Development economists don't believe SWID


# Geospatial GINI: Galimberti et al. (2020)
NL_GINI <- readxl::read_xlsx("data/GeospatialGinisData.xlsx", sheet = "Data") %>% qDT() %>% 
  fsubset(ISO == "ZAF", year = Year, wGini_L050) %>% fmutate(wGini_L050 = wGini_L050 * 100)
# Use SWIID after-tax income Gini-coefficients as a reference (Solt, 2016).

NL_GINI %>% ggplot(aes(x = year, y = wGini_L050)) + geom_line()

# Spatial Tax Panel
STP_MUN <- fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/Municipal_MedianIncome.csv") %>% fselect(-FTE) %>% 
  merge(fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/Municipal_Gini.csv"), by = .c(CAT_B, TaxYear)) %>% 
  fmutate(FTE = as.double(FTE))

STP_GINI <- STP_MUN %>% 
            fgroup_by(year = TaxYear) %>% 
            fsummarise(gini_between_mun = gini_wiki(MedianIncome)*100,
                       w_gini_between_mun = w_gini(MedianIncome, FTE)*100, 
                       gini_within_mun = fmean(gini)*100,
                       w_gini_within_mun = fmean(gini, FTE)*100, 
                       FTE = fsum(FTE), 
                       N = GRPN())            # fmutate(w_gini_adj = w_gini * (2 * FTE) / (N - 1))

STP_MUN %>% 
  fmutate(gini = gini * 100) %>% 
  ggplot(aes(x = TaxYear, y = gini, group = CAT_B, alpha = I(0.5))) + 
  geom_line() + guides(group = "none") +
  geom_line(aes(x = year, y = value, colour = Series), size = 1, alpha = I(1), inherit.aes = FALSE, 
            data = STP_GINI %>% gvr("year|within") %>% rm_stub("_within_mun", FALSE) %>% 
              melt("year", variable.name = "Series")) +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_bw(base_size = 14) + ylab("GINI") +
  ggtitle("STP3: GINI in 213 South African Municipalities")

dev.copy(pdf, "figures/STP3_SA_MUN_GINI_TimeSeries.pdf", width = 9.27, height = 5.83)
dev.off()

STP_GINI %>% gvr("year|gini") %>% melt("year") %>% 
  ggplot(aes(x = year, y = value, colour = variable)) + geom_line() +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  scale_color_brewer(palette = "Paired") +
  theme_bw(base_size = 14) + 
  ggtitle("STP3: GINI Between and Within 213 SA Municipalities")

dev.copy(pdf, "figures/STP3_SA_GINI_TimeSeries.pdf", width = 9.27, height = 5.83)
dev.off()


# All combined
SA_GINI_ALL <- SA_WID %>% fselect(series = variable, year, value) %>% mtt(series = paste("WID:", series)) %>% 
  rbind(SA_WB %>% fcompute(year = year(Date), series = "WB: SI_POV_GINI", value = SI_POV_GINI)) %>% 
  rbind(SA_SWID %>% melt("year", variable.name = "series") %>% mtt(series = paste("SWID:", series))) %>% 
  rbind(NL_GINI %>% fcompute(year = year, series = "DMSP-OLS: wGini_L050", value = wGini_L050)) %>% 
  rbind(STP_GINI %>% gvr("year|gini_within") %>% melt("year", variable.name = "series") %>% mtt(series = paste("STP3:", series)))
  

SA_GINI_ALL %>% 
  fsubset(year > 1990) %>% 
  frename(tools::toTitleCase) %>% 
  ggplot(aes(x = Year, y = Value, colour = Series)) + geom_line() + 
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  # scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) + 
  ggtitle("Different GINI Time Series for South Africa")

dev.copy(pdf, "figures/SA_GINI_TimeSeries.pdf", width = 9.27, height = 5.83)
dev.off()

SA_GINI_ALL %>% dcast(year ~ series) %>% pwcor()
  
#
### Optimization with Single kappa Objective ----------------------------
#

# swid_gini <- SA_SWID[between(year, 2014, 2022)]
# print(swid_gini)
np_pop_data_pos <- nl_pop_data %>% fsubset(pop > 0 & avg_rad > 0 & year == 2014, # year %in% swid_gini$year, 
                                           year, pop, avg_rad) # %>% fgroup_by(year)

nl_pop_data %>% 
  fgroup_by(nl_g0 = avg_rad > 0) %>% fselect(pop) %>% fsum()

pwcor(np_pop_data_pos)

objective <- function(k) {
    np_pop_data_pos %>% 
      fsummarise(nl_gini = w_gini(avg_rad^k, pop)*100) %$%
      # merge(swid_gini, by = "year") %$% 
      mean(abs(63 - nl_gini)) # World Bank GINI for 2014
}

result <- optimize(objective, c(0.01, 100))
# gini_disp: k = 1.334723, objective = 0.3401935
# gini_mkt: k = 1.822889, objective = 0.3002254

k <- 1.308966

sk_gini_res <- nl_pop_data %>% 
  fsubset(pop > 0 & avg_rad > 0) %>%
  fgroup_by(year) %>% 
  fsummarise(gini = gini_noss(avg_rad^k)*100, 
             w_gini = w_gini(avg_rad^k, pop)*100) 

sk_gini_res %$% ts(cbind(gini, w_gini), start = year[1]) %>% plot()

sk_gini_res %>% 
  fselect(-gini) %>% 
  merge(SA_WB %>% fcompute(year = year(Date), wb_gini = SI_POV_GINI), by = "year", all = TRUE) %>% 
  melt("year", na.rm = TRUE) %>% 
  frename(year = Year, value = Value, variable = Series) %>% 
  ggplot(aes(x = Year, y = Value, colour = Series)) + 
      geom_line() + 
      scale_y_continuous(limits = c(50,70), n.breaks = 10) +
      scale_x_continuous(n.breaks = 10) +
      scale_color_brewer(palette = "Set1") +
      theme_bw(base_size = 14) + 
      ggtitle("Calibrated (World Bank 2014) VIIRS-DNB Nightlights GINI Estimate")

dev.copy(pdf, "figures/WB_VIIRS_SA_GINI_TimeSeries.pdf", width = 9.27, height = 5.83)
dev.off()


#
### Computing Municipal GINIs ----------------------------
#

fastverse_extend(sf)
MunGeo <- st_read("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Shapefiles/MDB_Local_Municipal_Boundary_2018")
nl_pop_sf <- nl_pop_data %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)
system.time(ind <- st_contains(MunGeo, nl_pop_sf))
rm(nl_pop_sf); gc()

nl_pop_data$CAT_B <- NA_character_
for(i in seq_along(ind)) setv(nl_pop_data$CAT_B, ind[[i]], MunGeo$CAT_B[i], vind1 = TRUE)

wb_nl_gini_mun <- nl_pop_data %>% 
  fsubset(pop > 0 & avg_rad > 0 & !is.na(CAT_B)) %>% 
  fgroup_by(CAT_B, year) %>% 
  fsummarise(avg_rad = fmean(avg_rad),
             pop = fsum(pop),
             nl_gini = gini_noss(avg_rad^1.566813)*100, # Powers based on STP3 calibration
             nl_w_gini = w_gini(avg_rad^1.655078, pop)*100) 

STP_MUN_NL <- STP_MUN %>% rnm(TaxYear = year) %>% 
  merge(wb_nl_gini_mun, by = .c(CAT_B, year)) %>% 
  colorder(FTE, pop, avg_rad, pos = "after")

# populated_mun <- STP_MUN_NL %$% CAT_B[pop > 200000 & year == 2020]
populated_mun <- c("WC024", "CPT", "NMA", "WC044", "WC023", "BUF", "EC122", "EC121", 
"EC157", "EC155", "EC139", "EC153", "EC443", "KZN216", "EC441", 
"KZN225", "ETH", "KZN292", "KZN284", "MAN", "KZN237", "KZN282", 
"NC091", "KZN238", "FS194", "FS184", "KZN263", "KZN252", "NW403", 
"GT421", "MP307", "GT485", "NW405", "EKU", "JHB", "GT481", "MP312", 
"NW383", "MP313", "NW373", "MP315", "TSH", "NW372", "MP326", 
"NW371", "MP324", "LIM472", "MP316", "NW375", "LIM473", "LIM476", 
"MP325", "LIM355", "LIM354", "LIM333", "LIM367", "LIM332", "LIM331", 
"LIM344", "LIM345", "LIM343")

# Cross-sectional correlation
STP_MUN_NL %>% collap(~ CAT_B, na.rm = TRUE) %>% 
  fselect(med_inc = MedianIncome, FTE:nl_w_gini) %>% 
  pwcor() %>% print(digits = 3, show = "lower.tri") 

# Within-correlations
STP_MUN_NL %>% STD(~ CAT_B, na.rm = TRUE, stub = FALSE) %>% 
  fselect(med_inc = MedianIncome, FTE:nl_w_gini) %>% 
  pwcor() %>% print(digits = 3, show = "lower.tri") 

# All correlations: export
STP_MUN_NL %>% 
  # fsubset(CAT_B %in% populated_mun) %>%
  list(Overall = .,
       Between = collap(., ~ CAT_B, na.rm = TRUE), 
       Within = STD(., ~ CAT_B, na.rm = TRUE, stub = FALSE)) %>% 
  lapply(fselect, med_inc = MedianIncome, FTE:nl_w_gini) %>% 
  lapply(. %>% pwcor(P = TRUE) %>% print(digits = 3, return = TRUE, show = "lower.tri")) %>% 
  unlist2d("Trans", "Variable") %>% 
  roworderv(neworder = with(., radixorder(group(Variable)))) %>% 
  xtable::xtable() %>% print(booktabs = TRUE, include.rownames = FALSE)

STP_MUN_NL %>% 
  ggplot(aes(x = year, y = nl_w_gini, group = CAT_B, alpha = I(0.5))) + 
  geom_line() + guides(group = "none") +
  geom_line(aes(x = year, y = value, colour = Series), size = 1, alpha = I(1), inherit.aes = FALSE,
            data = STP_MUN_NL %>% gby(year) %>% smr(gini = fmean(nl_w_gini), w_gini = fmean(nl_w_gini, pop)) %>% 
              melt("year", variable.name = "Series")) +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_bw(base_size = 14) + ylab("GINI") +
  ggtitle("VIIRS-DNB: Pop Weighted GINI in 213 South African Municipalities")

dev.copy(pdf, "figures/WB_VIIRS_SA_MUN_GINI_TimeSeries.pdf", width = 9.27, height = 5.83)
dev.off()

#
### Optimization with Multiple Kappa Objective ----------------------------
#

# Different powers as in the paper
kappa <- c(0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0)

nl_pop_data_pos <- nl_pop_data %>% fsubset(pop > 0 & avg_rad > 0 & year %in% swid_gini$year) %>% fgroup_by(year) 
attr(nl_pop_data_pos, "kappa_mat") <- lapply(as.list(kappa), \(k) nl_pop_data_pos$avg_rad^k) %>% set_names(paste0("avg_rad_", kappa)) %>% qM()

objective <- function(kvec) {
  nl_pop_data_pos %>% 
    ftransform(avg_rad_w = attr(., "kappa_mat") %*% kvec) %>%
    fsummarise(nl_gini = w_gini(avg_rad_w, pop)*100) %>% qDT() %>%
    merge(swid_gini, by = "year") %$% 
    mean(abs(gini_disp - nl_gini))
}

mv_res <- optim(rep(1.1, length(kappa)), objective, 
             lower = rep(0.01, length(kappa)), 
             upper = rep(100, length(kappa)), 
             method = "L-BFGS-B")

  # -> No convergence !! (127 iterations, then: NORM OF PROJECTED GRADIENT <= PGTOL)
