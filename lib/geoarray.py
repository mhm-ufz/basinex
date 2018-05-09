#! /usr/bin/env python
# -*- coding: utf-8 -*-

# monkey patch GeoArray

import numpy as np
from ufz.geoarray import *


def setMask(grid, mask):
    data = np.array(grid, copy=True, subok=False)
    data[..., mask] = grid.fill_value
    return array(
        data=data,
        yorigin=grid.yorigin,
        xorigin=grid.xorigin,
        origin=grid.origin,
        cellsize=grid.cellsize,
        proj=grid.proj,
        fill_value=grid.fill_value,
        mode=grid.mode)


core.GeoArray.setMask = setMask
 
