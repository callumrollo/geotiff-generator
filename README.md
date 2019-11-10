# geotiff-generator

A python script to generate a geotiff image from bathymetry data. This geotiff can be uploaded to piloting softwre tools to give an idea of the depth of water.

### **Do not use this for navigational purposes!**

The author is a grad student not a seafarer, consult your local harbourmaster for charts, navigational advice and the like.

This script uses coarse, often out dated bathymetric datta and does **not** produce navigational charts, no accuracy is guaranteed or responsibility accepted for ay uses or misuses of this software.

### How to run the script:

Simply run the script with your Python interpreter of choice. I have included an anaconda environment file listing the necessary libraries.

The script defaults to interactive mode asking you what geographical area you wish to plot and where your bathy data is stored. You can specify most of the relevant arguments by passing them directly to **tiff_maker** 

You can supply your own lon, lat and bathy as numpy arrays or point the script to netcdf bathymetry files you have downloaded from [GEBCO](https://www.gebco.net/data_and_products/gridded_bathymetry_data/) or [EMODnet](https://portal.emodnet-bathymetry.eu/).

The script is very basic and be easily customised. Currently two colourpalettes are offered for the bathymap. I have not implemented a jet style colourmap. Please don't do it.
