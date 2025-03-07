# geotiff-generator

A python script to generate a geotiff image from bathymetry data. This geotiff can be uploaded to piloting software tools to give an idea of the depth of water.

**Test it out on Binder's cloud hosting**
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/callumrollo/geotiff-generator/master) may take a few minutes to load

When Binder has loaded, click on demo-try-me.ipynb

### **Do not use this for navigational purposes!**

The author is a PhD student not a seafarer. Consult your local harbourmaster for charts, navigational advice and the like.

This script uses coarse, often outdated bathymetric data and does **not** produce navigational charts. No accuracy is guaranteed or responsibility accepted for any uses or misuses of this software.

Please don't sue me, I'm poor.

### How to use geotiff generator:
1. Clone/download this repo with the green button
2. Install the packages listed in the environment.yml file (the easiset way is to [create an environment with Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)) 
3. Activate this environment
4. You can run the script `geotiff_gen.py` from a Python prompt or launch jupyter notebook and open the `demo-try-me.ipynb`

If you're unsure how to do either of these, you can check out the interactive demo with the Binder link above.

The script defaults to interactive mode asking you what geographical area you wish to plot and where your bathy data is stored. You can specify most of the relevant arguments by passing them directly to the function `tiff_maker`.

You can supply your own lon, lat and bathy as numpy arrays or point the script to netcdf bathymetry files you have downloaded from [GEBCO](https://www.gebco.net/data_and_products/gridded_bathymetry_data/) or [EMODnet](https://portal.emodnet-bathymetry.eu/).

The script is very basic and can be easily customised. Currently three colourpalettes are offered for the bathymap.

--------------

### Known Issues
#### These all stem from the installation of [gdal](https://gdal.org/) and how these tools are added to the PATH

##### 1. PROJ: proj_create_from_database: Cannot find proj.db

Python can't find the database of projections it needs to georefrence the tiff. This is a paths issue and indicates theat GDAL has not been installed properly

**Solutions:**

In linux I suggest appending anaconda3/envs/geotiff-generator/bin to your path in your shell rc file:

PATH=$PATH:/home/username/anaconda3/envs/geotiff-generator/bin

In Windows find the path (search for proj.db in your Anaconda3 folder) and append it to the path in your python console with:

setx PATH %PATH%;path-to-your-proj-folder

**Altering the PATH can have severe unintended consequencess** so take care. 

##### 2. ERROR 1: PROJ: proj_create_from_database: Open of /home/anaconda3/envs/geotiff-generator/share/proj failed

If using Pycharm the above error is sometimes generated. Not sure why this happens, running the file from the comand line or another IDE fixed it for me.
