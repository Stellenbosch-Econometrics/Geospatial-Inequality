## Geospatial (Remote Sensing) Inequality Estimates for South Africa

***Note:* This repo is still being cleaned up. **

This is exploratory work and summarized in the presentation. It depends on different data sources and several analysis scripts described below.

### Data

- [1.6km International Wealth Index Estimates](https://doi.org/10.7910/DVN/5OGWYM) by [Lee and Braithwaite (2022)](https://www.sciencedirect.com/science/article/pii/S0305750X22002182) and [2.4km Relative Wealth Index Extimates](https://dataforgood.facebook.com/dfg/tools/relative-wealth-index) by [Chi et al. (2022)](https://doi.org/10.1073/pnas.2113658119) for South Africa are included in the repo as csv files under `data/`.

- [463.8m VIIRS Stray Light Corrected Nighttime Day/Night Band Composites Version 1](https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_DNB_MONTHLY_V1_VCMSLCFG) from [Earth Observation Group](https://eogdata.mines.edu/products/vnl/) are available through [Google Earth Engine](https://earthengine.google.com/). This research uses annual median composites of the data. These can be created and downloaded, together with the monthly images, through `code/south_africa_nightlights.ipynb`. For the time being, the extracted [monthly](https://drive.google.com/drive/folders/1qjuSpBe2Xv2iqgoKzxc4crD2E2a94Ng0?usp=share_link) and [annual median composite](https://drive.google.com/drive/folders/18xI75APNFkUx4pcTfFdX8Orm36lcLzva?usp=share_link) images are available through my Google Drive. 

- [1km Population Estimates (UN-Adjusted)](https://hub.worldpop.org/geodata/listing?id=75) for South Africa are available from [WorldPop](https://hub.worldpop.org/). These are easily downloaded with the script `code/download_SA_pop.R`. 

- [Spatial Tax Panel v3.7](https://spatialtaxdata.org.za/). You can get it by filing a request to Nomonde Mathambo: nmathambo@hsrc.ac.za or contacting Dieter Fintel. It seems that data can also be downloaded [here](https://spatialtaxdata.org.za/download-data-filter-form). 

### Analysis Scripts

- `south_africa_inequality.R`: Produces spatial remotely sensed inequality estimates for South Africa inside 96km2 hexagons or 1km2 interpolations with 5 or 10km radius - based on IWI, RWI and Nightlights in 2020. 

- `NA_NL_GINI_TS.R`: Produces time series estimates of inequality in South Africa based on annual nightlights median composite images since 2014. 

