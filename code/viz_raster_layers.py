# %% 
import os
import numpy as np
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


# %%
