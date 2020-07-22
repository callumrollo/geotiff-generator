from geotiff_gen import gebco_subset, emod_subset, argnearest
import math
import numpy as np

path_to_gebco = '/media/callum/storage/Documents/global_datasets/GEBCO/'
path_to_emod = '/media/callum/storage/Documents/global_datasets/emod_netcdf/'


def test_gebco_subset():
    lon_sub, lat_sub, bathy_sub = gebco_subset(path_to_gebco, [0, 1, -10, -9], True)
    assert math.isclose(lon_sub[0], -10, rel_tol=1e-2)
    assert math.isclose(lat_sub[-1], 1, rel_tol=1e-2)
    assert math.isclose(np.nanmean(bathy_sub), -4800, rel_tol=1e-1)


def test_emod_subset():
    # from 1 tile
    __ = emod_subset(path_to_emod, [53, 54, 2, 3])
    # from 2 zonal tiles
    __ = emod_subset(path_to_emod, [53, 54, 2, 4])
    # from 2 meridional tiles
    __ = emod_subset(path_to_emod, [52, 54, 2, 3])
    # from 4 tiles
    lon_sub, lat_sub, bathy_sub = emod_subset(path_to_emod, [52, 54, 2, 4])
    assert lon_sub[-1] - lon_sub[0] > 2.0
    assert lat_sub[-1] - lat_sub[0] > 2.0
    assert math.isclose(np.nanmean(bathy_sub), -34, rel_tol=1e-1)


def test_argnearest():
    points = np.arange(1, 10)
    assert argnearest(points, 3.3) == 2
