from geotiff_gen import argnearest

import numpy as np

def test_argnearest():
    points = np.arange(1,10)
    assert argnearest(points,3.3) == 2