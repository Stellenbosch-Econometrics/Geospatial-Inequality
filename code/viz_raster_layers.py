# %% 
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal
os.getcwd()

# %% 
minlon = 16; minlat = -35; maxlon = 33; maxlat = -22
bbox = (minlon, maxlat, maxlon, minlat)

# %%
data_path = "../data/south_africa_viirs_dnb_nightlights_v1_vcmslcfg_annual_median_composite"
images = sorted(os.listdir(data_path))
images

# %% 
fig, axs = plt.subplots(ncols=3, nrows=3, figsize=(10, 10))
fig.suptitle('VIIRS DNB Nightlights V1', fontsize=16)
fig.tight_layout()
# compress the vertical space between subplots
fig.subplots_adjust(top=0.75, hspace=0, wspace=0)

for i in range(len(images)):
    # Using GDAL virtual file system: https://gis.stackexchange.com/questions/264793/crop-raster-in-memory-with-python-gdal-bindings
    gdal.Translate('/vsimem/clip.tif', os.path.join(data_path, images[i]), projWin=bbox)
    img = gdal.Open('/vsimem/clip.tif')
    img_arr = np.log10(img.GetRasterBand(1).ReadAsArray())
    # replace infinite values with nan
    img_arr[img_arr == -np.inf] = np.nan
    del img
    # axs[int(i/3), i%3].set_aspect('equal')
    axs[int(i/3), i%3].margins(0) 
    axs[int(i/3), i%3].axis('off')
    axs[int(i/3), i%3].imshow(img_arr, cmap= "magma", interpolation='nearest')
    axs[int(i/3), i%3].set_title(images[i][:4] + " [{}, {}]".format(int(10**np.nanmin(img_arr)), int(10**np.nanmax(img_arr))))

plt.tight_layout()
plt.savefig('../figures/nightlights_layers.png', dpi = 200)
plt.show()

# %%
data_path = "../data/WPOP_SA_1km_UNadj"
images = sorted(os.listdir(data_path))
images

# %% 
fig, axs = plt.subplots(ncols=3, nrows=3, figsize=(10, 10))
fig.suptitle('WorldPop Population (UN Adjusted)', fontsize=16)
fig.tight_layout()
# compress the vertical space between subplots
# fig.subplots_adjust(top=0.75, hspace=0, wspace=0)

for i in range(9):
    if i < len(images):
        # Using GDAL virtual file system: https://gis.stackexchange.com/questions/264793/crop-raster-in-memory-with-python-gdal-bindings
        gdal.Translate('/vsimem/clip.tif', os.path.join(data_path, images[i]), projWin=bbox)
        img = gdal.Open('/vsimem/clip.tif')
        img_arr = np.log10(img.GetRasterBand(1).ReadAsArray())
        # replace infinite values with nan
        img_arr[img_arr == -np.inf] = np.nan
        del img
    # axs[int(i/3), i%3].set_aspect('equal')
    axs[int(i/3), i%3].margins(0) 
    axs[int(i/3), i%3].axis('off')
    if i < len(images):
        axs[int(i/3), i%3].imshow(img_arr, cmap= "viridis", interpolation='nearest')
        axs[int(i/3), i%3].set_title(images[i][8:12] + " [{}, {}]".format(int(10**np.nanmin(img_arr)), int(10**np.nanmax(img_arr))))

plt.tight_layout()
plt.savefig('../figures/population_layers.png', dpi = 200)
plt.show()


# %% Now we will plot the estimated wealth index and the RWI
SA_IWI = pd.read_csv("../data/SA_IWI.csv").rename({"estimated_IWI" : "IWI"}, axis=1)
SA_IWI = SA_IWI[SA_IWI.lon.between(minlon, maxlon) & SA_IWI.lat.between(minlat, maxlat)]
SA_RWI = pd.read_csv("../data/SA_RWI.csv").rename({"longitude" : "lon", "latitude" : "lat", "rwi" : "RWI"}, axis=1)
SA_RWI = SA_RWI[SA_RWI.lon.between(minlon, maxlon) & SA_RWI.lat.between(minlat, maxlat)]

# %%
import geopandas as gpd
south_africa = gpd.read_file("../data/spatial_tax_panel/Spatial_Tax_Panel_v3/Shapefiles/MDB_Local_Municipal_Boundary_2018/MDB_Local_Municipal_Boundary_2018.shp") \
                  .clip_by_rect(xmin = minlon, ymin = minlat, xmax = maxlon, ymax = maxlat) 
# %% 
# IWI
fig, ax = plt.subplots(figsize = (12, 12))
# fig.suptitle('International Wealth Index', fontsize=16)
ax.set_aspect('equal')
ax.margins(0) 
ax.set_axis_off() # Only Bbox: ax.set_frame_on(False)
south_africa.plot(ax = ax, linewidth = 0.05, color='grey', edgecolor='black') # black
plt.scatter(x = SA_IWI.lon, y = SA_IWI.lat, c = SA_IWI.IWI, marker=".", 
            s = 3, cmap="turbo", edgecolors='none', vmin=0, vmax=100)
cb = plt.colorbar(shrink = 0.5, pad = 0.02)
cb.outline.set_color(None)
cb.ax.set_title('IWI')
plt.savefig("../figures/SA_IWI.pdf", 
            dpi = 300, format="pdf", bbox_inches="tight")
plt.savefig("../figures/SA_IWI.png", 
            dpi = 300, format="png", bbox_inches="tight")
plt.show()
# %% 
# RWI
fig, ax = plt.subplots(figsize = (12, 12))
ax.set_aspect('equal')
ax.margins(0) 
ax.set_axis_off() # Only Bbox: ax.set_frame_on(False)
south_africa.plot(ax = ax, linewidth = 0.05, color='grey', edgecolor='black') # black
plt.scatter(x = SA_RWI.lon, y = SA_RWI.lat, c = SA_RWI.RWI, marker="s", 
            s = 0.57, cmap="turbo", edgecolors='none') # , vmin=0, vmax=100
cb = plt.colorbar(shrink = 0.5, pad = 0.02)
cb.outline.set_color(None)
cb.ax.set_title('RWI')
plt.savefig("../figures/SA_RWI.pdf", 
            dpi = 300, format="pdf", bbox_inches="tight")
plt.savefig("../figures/SA_RWI.png", 
            dpi = 300, format="png", bbox_inches="tight")
plt.show()
# %%
## NL20
data_path = "../data/south_africa_viirs_dnb_nightlights_v1_vcmslcfg_annual_median_composite"
gdal.Translate('/vsimem/clip.tif', os.path.join(data_path,  '2020.tif'), projWin=bbox)
img = gdal.Open('/vsimem/clip.tif')
img_arr = np.log10(img.GetRasterBand(1).ReadAsArray())
# replace infinite values with nan
img_arr[img_arr == -np.inf] = np.nan
del img
# %%
fig, ax = plt.subplots(figsize = (12, 12))
ax.set_aspect('equal')
ax.margins(0) 
ax.set_axis_off() # Only Bbox: ax.set_frame_on(False)
plt.imshow(img_arr, cmap="turbo", interpolation='nearest')
# south_africa.plot(ax = ax, linewidth = 0.05, color='none', edgecolor='black') # black
cb = plt.colorbar(shrink = 0.5, pad = 0.02)
cb.outline.set_color(None)
cb.ax.set_title('NL20')
plt.savefig("../figures/SA_NL20.png", 
            dpi = 200, format="png", bbox_inches="tight")
plt.show()
# %%
## POP20
data_path = "../data/WPOP_SA_1km_UNadj"
gdal.Translate('/vsimem/clip.tif', os.path.join(data_path,  'zaf_ppp_2020_1km_Aggregated_UNadj.tif'), projWin=bbox)
img = gdal.Open('/vsimem/clip.tif')
img_arr = np.log10(img.GetRasterBand(1).ReadAsArray())
# replace infinite values with nan
img_arr[img_arr == -np.inf] = np.nan
del img
# %%
fig, ax = plt.subplots(figsize = (12, 12))
ax.set_aspect('equal')
ax.margins(0) 
ax.set_axis_off() # Only Bbox: ax.set_frame_on(False)
plt.imshow(img_arr, cmap="turbo", interpolation='nearest')
# south_africa.plot(ax = ax, linewidth = 0.05, color='none', edgecolor='black') # black
cb = plt.colorbar(shrink = 0.5, pad = 0.02)
cb.outline.set_color(None)
cb.ax.set_title('POP20')
plt.savefig("../figures/SA_POP20.png", 
            dpi = 200, format="png", bbox_inches="tight")
plt.show()
# %%
