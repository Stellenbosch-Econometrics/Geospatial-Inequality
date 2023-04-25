## Geospatial (Remote Sensing) Inequality Estimates for South Africa

This is exploratory work and summarized in the presentation. It depends on different data sources and several analysis scripts described below.

### Data

- [1.6km International Wealth Index Estimates](https://doi.org/10.7910/DVN/5OGWYM) by [Lee and Braithwaite (2022)](https://www.sciencedirect.com/science/article/pii/S0305750X22002182) and [2.4km Relative Wealth Index Extimates](https://dataforgood.facebook.com/dfg/tools/relative-wealth-index) by [Chi et al. (2022)](https://doi.org/10.1073/pnas.2113658119) for South Africa are included in the repo as csv files under `data/`.

- [463.8m VIIRS Stray Light Corrected Nighttime Day/Night Band Composites Version 1](https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_DNB_MONTHLY_V1_VCMSLCFG) from [Earth Observation Group](https://eogdata.mines.edu/products/vnl/) are available through [Google Earth Engine](https://earthengine.google.com/). This research uses annual median composites of the data. These can be created and downloaded, together with the monthly images, through `code/south_africa_nightlights.ipynb`. For the time being, the extracted [monthly](https://drive.google.com/drive/folders/1qjuSpBe2Xv2iqgoKzxc4crD2E2a94Ng0?usp=share_link) and [annual median composite](https://drive.google.com/drive/folders/18xI75APNFkUx4pcTfFdX8Orm36lcLzva?usp=share_link) images are available through my Google Drive. 

- An improved [Annual Median Composite (V2.1)](https://eogdata.mines.edu/products/vnl/#annual_v2) of the VIIRS Nightlights, also at 463.8m resolution, for South Africa is included under `data/`. The corresponding global image was downloaded from [here](https://eogdata.mines.edu/nighttime_light/annual/v21/2021/) and cropped to South Africa using a rectangular bounding box. This image is used for all geospatial analysis, V1 images are only used to generate time series. A variant of this image with zeros replaced by `NA`'s is included for visualization purposes. 

- [1km Population Estimates (UN-Adjusted)](https://hub.worldpop.org/geodata/listing?id=75) for South Africa are available from [WorldPop](https://hub.worldpop.org/). These are easily downloaded with the script `code/download_SA_pop.R`. 

- Various Official GINI Estimates: `data/SA_GINI` has estimates from the [World Inequality Databse](https://wid.world/) (extracted using the STATA API) and [Standardized World Inequality Database](https://fsolt.org/swiid/). Furthermore `data/GeospatialGinisData.xlsx` has nightlights-based GINI estimates for all countries from 1992-2013, based on the old DMSP-OLS satellite and the work of [Galimberti et al. (2020)](https://www.aut.ac.nz/__data/assets/pdf_file/0010/399394/working-paper-20_07.pdf), downloaded from [here](https://www.ciesin.columbia.edu/data/global-geospatial-inequality/). 

- [Spatial Tax Panel v3.7](https://spatialtaxdata.org.za/). You can get it by filing a request to Nomonde Mathambo: nmathambo@hsrc.ac.za or contacting Dieter Fintel. It seems that data can also be downloaded [here](https://spatialtaxdata.org.za/download-data-filter-form). Results including municipal and hexagonal shapes are exported - together with the relevant STP3 data, as [GeoPackage](https://www.geopackage.org/) databases under `results/` (thus it is actually not really necessary to get the STP3, unless you want to replicate the results on a newer version). 

### Analysis Scripts

- `SA_NL_GINI_TS.R`: Produces time series estimates of inequality in South Africa based on annual nightlights median composite images since 2014. It also computes and explores municipal nightlights time series, which are however deemed useless.

- `SA_inequality.R`: Produces spatial remotely sensed inequality estimates for South Africa inside 96km2 hexagons or 1km2 interpolations with 5 or 10km radius - based on IWI, RWI and Nightlights in 2020. 

- `explore_SA_inequality.R`: Explores the inequality estimates computed in `SA_inequality.R` using correlations and graphs, and relates them to the Uber Hexagons of the Spatial Tax Panel v3.7. 

- `spatial_tax_panel.R`: Explores the Spatial Tax Panel v3.7 and joins it with the RWI, IWI, Nightlights and Population, from which alternative municipal GINI estimates are produces. The estimates are compared using correlations and graphs. 

- `explore_spatial_tax_panel.R`: Explores the inequality estimates computed in `spatial_tax_panel.R` using correlations and graphs. 

- `viz_raster_layers.py`: Python script to plot raster layers and nightlights/wealth estimates at high resolution using Matplotlib. 

### Results

- The `results/` folder contains all the GINI estimates produced by various scripts. The result involving 96km2 hexagons is saved as a [GeoPackage](https://www.geopackage.org/) database, which contains the geometry and can be read from many softwares. The folder also contains a [QGIS](https://qgis.org/en/site/) project to visualize these different estimates. 

- The `figures/` folder contains various graphs and figures, many of which are included in the presentation. 

- The `presentation/` folder contains a Beamer presentation of the results, delivered at the BBL seminar in the Stellenbosch Economics Department on April 18, 2023. The seminar was recorded. 

### Further Notes

- Every script should be evaluated on a fresh R session, in particular settings that optimize some of the libraries used such as `set_collapse(nthreads = 4, na.rm = FALSE, sort = FALSE)` are not to be used on all scripts. 
