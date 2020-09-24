#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:35:59 2019

@author: callum

An interactive script for making geotiff images from bathymetry data
"""
import os
import sys
from tkinter import Tk
from tkinter.filedialog import askopenfilename, askdirectory
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import copy
from pathlib import Path
from osgeo import gdal, osr



def main():
    tiff_maker()


def tiff_maker(filename='', lon=[], lat=[], bathy=[], extent=[], bathy_type='', bathy_folder_path='', theme='', min_depth='usr', bathy_nc=False):
    """
    Interactive function for generating geotiffs from one of three sources. Can be called empty for
    interactive use or you can specify all the arguments you will nedd with the kwargs as follows

    :param filename: filename for the geotiff (do not specify an extension) always required
    :param lon: numpy array or list of lon (if using own data)
    :param lat: numpy array or list of lat (if using own data)
    :param bathy: numpy array or 2D list. dimensions must match lon and lat (if using own data)
    :param extent: list with four items which are extent of desired geotiff [South, North, West, East]
    e.g. [49. 50.5, -5, 2] (if using gebco or emodnet)
    :param bathy_type: 'e' for emodnet, 'g' for GEBCO
    :param bathy_folder_path: path to the folder with your bathymetry netCDFs (if using GEBCO or emodnet)
    :param theme: colour theme of the geotiff 'y' for dark 'n' for light 'g' for greyscale
    :param min_depth: the minimum depth, all water shallower than this coloured red
    :param bathy_nc set to True to write a netcdf file fo the subset bathymetry (currently supports only GEBCO)
    """
    if not filename:
        filename = input("Enter a name for your bathymetry geotiff ")
    if len(lon) != 0:
        bathy_to_tiff(lon, lat, bathy, filename, theme, min_depth)
        return
    if not extent:
        extent = []
        extent.append(float(input("Enter the southern limit of your desired bathy ")))
        extent.append(float(input("Enter the northern limit of your desired bathy ")))
        extent.append(float(input("Enter the western limit of your desired bathy ")))
        extent.append(float(input("Enter the eastern limit of your desired bathy ")))
    if extent[1]<extent[0]:
        print(f"northern limit ({extent[1]}) must be greater than southern limit ({extent[0]}). Aborting")
        sys.exit()
    if extent[3]<extent[2]:
        print(f"eastern limit ({extent[3]}) must be greater than western limit ({extent[2]}). Aborting")
        sys.exit()
    if bathy_type in ['e', 'g']:
        bathy_selec=bathy_type
    else:
        bathy_selec = input("What bathy are you using? [g]ebco, [e]modnet or [o]ther ")
    if bathy_selec.lower() == 'o':
        print("call  tiff_maker with your lon, lat and bathy arrays\n"
              "e.g. tiff_generator(filename='my_cool_tiff_pic', lon=<lon_array>, lat=<lat_array>, bathy=<bathy_array>)")
        return
    if bathy_selec.lower() == 'g':
        if bathy_folder_path:
            gebco_path = bathy_folder_path
        else:
            selection = input(r"Enter the path to a GEBCO netcdf file. Press 'f' for a file selector app ")
            if selection.lower() == "f":
                print("take me to your file")
                Tk().withdraw() 
                gebco_path = askopenfilename()
            else:  
                gebco_path = selection
        print("Using GEBCO bathy from " + gebco_path)
        lon, lat, bathy = gebco_subset(gebco_path, extent, bathy_nc=bathy_nc)
        bathy_to_tiff(lon, lat, bathy, filename, theme, min_depth)
        return
    if bathy_selec.lower() == 'e':
        if bathy_folder_path:
            emod_path = bathy_folder_path
        else:
            selection = input(r"Enter the path to the folder with your emodnet file(s) (E3.dtm etc). Press 'f' for a folder selector app  ")
            if selection.lower() == "f":
                print("take me to your folder")
                Tk().withdraw() 
                emod_path = askdirectory()
            else:  
                emod_path = selection
        lon, lat, bathy = emod_subset(emod_path, extent)
        bathy_to_tiff(lon, lat, bathy, filename, theme, min_depth)
        return


def bathy_to_tiff(lon_vals, lat_vals, bathy_vals, filename, theme, min_depth):
    """
    The script that makes the geotiff from numpy arrays of lon, lat and bathy
    :param lon_vals: 1D np.array of lon values
    :param lat_vals: 1D np.array of lat values
    :param bathy_vals: 2D np.array of bathy data, must conform with lon and lat
    :param filename: base string for the filename, do not specify an extension
    """
    # Set all nans to height 0.0 m
    bathy_vals[np.isnan(bathy_vals)] = 0.0

    # set all land to height 0.0 m
    bathy_vals[bathy_vals > 0.0] = 0.0

    # If needed, flip the y axis for writing the file to tif (grid must be in shape of map so starting from the NW corner)
    if lat_vals[-1] > lat_vals[0]:
        lat_vals = lat_vals[::-1]
        bathy_vals = bathy_vals[::-1, :]

    #  Initialize the image
    image_size = (len(lat_vals), len(lon_vals))

    #  Create red channel, starts set to 0
    r_pixels = np.zeros(image_size, dtype=np.uint8)

    # Gradient in blue and green  channels from bathy
    g_pixels = 255 * (1 - np.abs(bathy_vals / np.nanmin(bathy_vals)))
    b_pixels = 255 * (1 - np.abs(bathy_vals / np.nanmin(bathy_vals)))

    if not theme:
        print("Would you like the geotiff in a dark theme?\n"
              "[y]es please, I love dark themed UI elements\n"
              "[n]ah just melt my eyeballs out\n"
              "[g]reyscale for true bathymetry")
        while True:
            theme = input("")
            if theme.lower() not in ('y', 'n', 'g'):
                print("\nPlease enter y, n or g\n")
            else:
                break

    theme_color = [0, 80, 0, 120]
    if theme == 'n':
        theme_color = [255, 255, 255, 255]
    # Set land to theme color
    r_pixels[bathy_vals >= 0] = theme_color[0]
    g_pixels[bathy_vals >= 0] = theme_color[1]
    b_pixels[bathy_vals >= 0] = theme_color[2]

    # Set all heights on a greyscale for true scaled bathymetry and topography
    if theme == 'g':
        r_pixels = g_pixels = b_pixels = 255 * (np.abs(bathy_vals / np.nanmin(bathy_vals)))
    else:
        if type(min_depth) == str:
            print("What depth (m) would you like the shallow warning red set? (0 for no shallow warning) ")
            min_depth = float(input(""))

        # Set bathy shallower than user desired value to red (intensity depends on theme chosen)
        r_pixels[np.logical_and(bathy_vals > -np.abs(min_depth), bathy_vals < 0)] = theme_color[3]
        g_pixels[np.logical_and(bathy_vals > -np.abs(min_depth), bathy_vals < 0)] = 0
        b_pixels[np.logical_and(bathy_vals > -np.abs(min_depth), bathy_vals < 0)] = 0

    print("Generating geotiff...")
    # set geotransform
    nx = image_size[1]
    ny = image_size[0]
    xmin, ymin, xmax, ymax = [min(lon_vals), min(lat_vals), max(lon_vals), max(lat_vals)]
    xres = (xmax - xmin) / float(nx)
    yres = (ymax - ymin) / float(ny)
    geotransform = (xmin, xres, 0, ymax, 0, -yres)

    # create the 3-band raster file

    dst_ds = gdal.GetDriverByName('GTiff').Create(filename + '.tif', nx, ny, 3, gdal.GDT_Byte)

    srs = osr.SpatialReference()  # establish encoding
    srs.ImportFromEPSG(4326)  # Import the WGS84 datum
    dst_ds.SetGeoTransform(geotransform)  # specify coords
    dst_ds.SetProjection(srs.ExportToWkt())  # export coords to file
    dst_ds.GetRasterBand(1).WriteArray(r_pixels)  # write r-band to the raster
    dst_ds.GetRasterBand(2).WriteArray(g_pixels)  # write g-band to the raster
    dst_ds.GetRasterBand(3).WriteArray(b_pixels)  # write b-band to the raster
    dst_ds.FlushCache()  # write to disk
    dst_ds = None  # clean up
    print('Made geotiff file at: ' + str(Path(os.getcwd()) / (filename + '.tif')))


def argnearest(items, pivot):
    near_item = min(items, key=lambda x: abs(x - pivot))
    for i in range(len(items)):
        if items[i] == near_item:
            return i


def gebco_subset(path_to_folder_str, extent, bathy_nc):
    """
    Extracts bathy data from a global GEBCO .nc file from an area specified by the use
    :param path_to_folder: string of path to the folder, specified by user
    :param extent: list with four items which are extent of desired geotiff [South, North, West, East]
    e.g. [49. 50.5, -5, 2] (if using gebco or emodnet)
    :return: numpy arrays of lon, lat and bathymetry
    """
    print('Fetching GEBCO data...')
    path_to_folder = Path(path_to_folder_str)
    if path_to_folder.is_file():
        gebco = xr.open_dataset(path_to_folder)
    else:
        path_to_gebco = list(Path(path_to_folder).joinpath().glob("*.nc"))
        if not path_to_gebco:
            print('No netcdf files found in location supplied. Check that you pointed to a .nc file or a folder containing one. Aborting')
            exit(1)
        gebco = xr.open_dataset(path_to_gebco[0])
    print("Subsettting GEBCO data")
    subset = gebco.sel(lon=slice(extent[2], extent[3]), lat=slice(extent[0], extent[1]))
    "print GEBCO bathy fetch successful"
    if bathy_nc==True:
        ## To save our bathymetry data
        subset.to_netcdf(Path(os.getcwd())/'bathy_subset.nc')
        print('bathy subset written at ' + str(Path(os.getcwd())/'bathy_subset.nc'))
    return np.array(subset.lon), np.array(subset.lat), np.array(subset.elevation)


def emod_subset(path_to_files, extent):
    """
    Selects the EMODnet bathymetry tiles that cover the user defined area and
    combines them for a seamless bathymetry. Returns lon, lat and bathy.
    :param path_to_files: string of path to the folder containing EMODnet .dtm files, specified by user
    :param extent: list with four items which are extent of desired geotiff [South, North, West, East]
    e.g. [49. 50.5, -5, 2] (if using gebco or emodnet)
    :return: numpy arrays of lon, lat and bathymetry
    """
    print("Searching recursivly for emod *.dtm bathy files in " + path_to_files)
    tiles = Path(path_to_files).joinpath().glob("**/*.dtm")
    tiles_check = Path(path_to_files).joinpath().glob("**/*.dtm")
    if not list(tiles_check):
        print(
            'No netcdf files found in supplied folder. Check that it is a complete folder path (not a file). Aborting')
        sys.exit()
    unorder_list = []
    for item in tiles:
        unorder_list.append(str(item))
    tile_list = np.sort(unorder_list)
    tile_name = copy.deepcopy(tile_list)
    extents = np.empty((len(tile_list), 4))

    for tile in range(len(tile_list)):
        tile_name[tile] = tile_list[tile][-11:-9]
        ds = xr.open_dataset(tile_list[tile])
        attributes = ds.attrs
        # Get the locations of the bottom left and top right corners
        extents[tile, 2] = attributes["Longitude_BL"]
        extents[tile, 3] = attributes["Longitude_TR"]
        extents[tile, 0] = attributes["Latitude_BL"]
        extents[tile, 1] = attributes["Latitude_TR"]

    # Put extents into a pretty dataframe
    tile_extents = pd.DataFrame(
        data=extents, columns=["South", "North", "West", "East"], index=tile_name
    )
    S_in, N_in, W_in, E_in = extent
    # Add a 0.1 degree buffer around selected area
    S = S_in - 0.1
    N = N_in + 0.1
    W = W_in - 0.1
    E = E_in + 0.1
    relevant_tiles = []
    vertices = [[W, S], [E, S], [W, N], [E, N]]
    # Find tiles within user defined region
    for tile in tile_name:
        this_tile = tile_extents.loc[tile]
        for vertex in range(4):
            if np.logical_and(
                    vertices[vertex][0] > this_tile.loc["West"],
                    vertices[vertex][0] < this_tile.loc["East"],
            ):
                if np.logical_and(
                        vertices[vertex][1] > this_tile.loc["South"],
                        vertices[vertex][1] < this_tile.loc["North"],
                ):
                    relevant_tiles.append(tile)
                    break
    # Sort the tiles by column first ready for patching
    relevant_tiles.sort(key=lambda x: x[-1])
    tile_paths = []
    for tile in relevant_tiles:
        tile_paths.append(tile_list[tile_name == tile][0])
    if len(relevant_tiles)==0:
        print("Area specified has no overlap with EMODnet bathy supplied. Aborting")
        sys.exit()
    print(f"Bathymetry data contained in {len(relevant_tiles)} files. Fetching...")

    # Extracting the relevant bathymetry data from the first tile
    tile = xr.open_dataset(tile_paths[0])
    lon_tile = tile.COLUMNS[np.logical_and(tile.COLUMNS > W, tile.COLUMNS < E)]
    lat_tile = tile.LINES[np.logical_and(tile.LINES > S, tile.LINES < N)]
    bathy_tile = tile.DEPTH.sel(COLUMNS=lon_tile, LINES=lat_tile)
    # Now to add bathy data from more tiles if needed

    if len(relevant_tiles) > 1:
        base_tile = relevant_tiles[0]
        print(f"base tile {base_tile}")
        for i, next_tile in enumerate(relevant_tiles[1:]):
            print(f"adding tile {next_tile}")
            if next_tile[1] == base_tile[1]:
                # For tiles in one column
                print("patching NS...")
                add_tile = xr.open_dataset(tile_paths[i + 1])
                add_lat = add_tile.LINES[
                    np.logical_and(add_tile.LINES > S, add_tile.LINES < N)
                ]
                add_lat_no_overlap = add_lat[add_lat < np.nanmin(lat_tile)]
                if len(add_lat_no_overlap) == 0:
                    continue
                add_bathy = add_tile.DEPTH.sel(
                    COLUMNS=lon_tile, LINES=add_lat_no_overlap
                )
                lat_tile = np.concatenate((add_lat_no_overlap, lat_tile))
                bathy_tile = np.concatenate((add_bathy, bathy_tile))

            elif next_tile[0] == base_tile[0]:
                # If tiles are from more than one column, initiates a new column
                print("start new column...")
                column_base_tile = xr.open_dataset(tile_paths[i + 1])
                column_lon_tile = column_base_tile.COLUMNS[
                    np.logical_and(column_base_tile.COLUMNS > np.nanmax(lon_tile), column_base_tile.COLUMNS < E)]
                column_lat_tile = column_base_tile.LINES[
                    np.logical_and(column_base_tile.LINES > S, column_base_tile.LINES < N)]
                column_bathy_tile = column_base_tile.DEPTH.sel(COLUMNS=column_lon_tile, LINES=column_lat_tile)
            else:
                # Adds to the new column
                print("patching NS...")
                add_tile = xr.open_dataset(tile_paths[i + 1])
                add_lat = add_tile.LINES[
                    np.logical_and(add_tile.LINES > S, add_tile.LINES < N)
                ]
                add_lat_no_overlap = add_lat[add_lat < np.nanmin(column_lat_tile)]
                add_bathy = add_tile.DEPTH.sel(
                    COLUMNS=column_lon_tile, LINES=add_lat_no_overlap
                )
                column_lat_tile = np.concatenate((add_lat_no_overlap, column_lat_tile))
                column_bathy_tile = np.concatenate((add_bathy, column_bathy_tile))
            if not 'column_lat_tile' in locals():
                # Exit if we haven't started a second column yet( happens with zero overlap row)
                continue
            if np.logical_and(next_tile[1] != base_tile[1], len(lat_tile) == len(column_lat_tile)):
                # Combines the columns
                print("combining columns")
                add_lon_no_overlap = column_lon_tile[column_lon_tile > np.nanmax(lon_tile)]
                lon_tile = np.concatenate((lon_tile, add_lon_no_overlap))
                bathy_tile = np.concatenate((bathy_tile, column_bathy_tile), axis=1)
    if type(lon_tile) != np.ndarray:
        lon_tile = lon_tile.values
    if type(lat_tile) != np.ndarray:
        lat_tile = lat_tile.values
    if type(bathy_tile) != np.ndarray:
        bathy_tile = bathy_tile.values
    print('Bathy parsed successfully')
    return lon_tile, lat_tile, bathy_tile


if __name__ == '__main__':
    main()

