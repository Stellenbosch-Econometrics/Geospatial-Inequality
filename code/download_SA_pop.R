########################################################
# Download WorldPop UN-Adjusted 1km2 Population for SA
# Source: https://hub.worldpop.org/geodata/listing?id=75
########################################################

pop_year <- function(y) paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", y, "/ZAF/zaf_ppp_", y, "_1km_Aggregated_UNadj.tif")

for (y in 2014:2020) {
  download.file(pop_year(y), mode = "wb", 
                destfile = paste0("data/WPOP_SA_1km_UNadj/zaf_ppp_", y, "_1km_Aggregated_UNadj.tif"))
}
