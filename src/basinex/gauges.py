# -*- coding: utf-8 -*-

import csv

import numpy as np


class Gauge(object):
    def __init__(
        self, id, y=None, x=None, size=None, path=None, varname=None, lat_fix=False
    ):
        self.id = id
        self.y = y
        self.x = x
        self.size = float(size) / np.cos(np.deg2rad(float(y))) if lat_fix else size
        self.path = path
        self.varname = varname

    def __str__(self):
        return str(
            {
                "id": self.id,
                "y": self.y,
                "x": self.x,
                "size": self.size,
                "path": self.path,
                "varname": self.varname,
            }
        )


def readGauges(fname, lat_fix=False):
    keysets = {
        tuple(sorted(["id", "size", "y", "x"])),
        tuple(sorted(["id", "path", "varname"])),
    }
    with open(fname) as f:
        reader = csv.DictReader(f, delimiter=";")
        out = []
        for gauge in reader:
            if tuple(sorted(gauge.keys())) not in keysets:
                raise RuntimeError(
                    "Header in gauges look up table must be in {:}".format(
                        list(keysets)
                    )
                )
            out.append(Gauge(lat_fix=lat_fix, **gauge))
    return out


def matchFlowacc(gauge, facc, max_distance, max_error, scaling_factor=1):
    """
    Input:
        GeoGrid or sibling, Numeric, Numeric
    Purpose:
        - Moves the Instance onto the nearest cell in the input
          flowaccumulation map (facc) based on the catchment_area.
        - The search will be done within the rectangele
          coordinate - max_distance to coordinate + max_distance
        - The method minimizes the error between the attribute
          catchment_area and the number of cells acuumulated in the
          grid. If the resulting difference is larger than
          max_error False will be returned.
        - The argument scaling factor allows to handle different
          horizontal resolutions of the information given. A factor
          to convert from map units to km^2 of the catchment area
    """
    # print "gauge before:", gauge.y, gauge.x
    y = float(gauge.y)
    x = float(gauge.x)
    size = float(gauge.size)
    bbox = {
        "ymin": y - max_distance,
        "ymax": y + max_distance,
        "xmin": x - max_distance,
        "xmax": x + max_distance,
    }
    grid = facc.shrink(**bbox).astype("float32")
    grid *= abs(
        (grid.cellsize[0] * scaling_factor) * (grid.cellsize[1] * scaling_factor)
    )

    gauge_y_idx, gauge_x_idx = grid.indexOf(y, x)

    # find all possible river cells in the nearest possible distance
    error = 0
    river_cells_y = []
    while len(river_cells_y) == 0 and error < max_error:
        river_cells_y, river_cells_x = np.where(
            (grid.data >= size * (1 - error)) & (grid.data <= size * (1 + error))
        )
        error += 0.01

    if error <= max_error:
        # the closest river cell
        nn = (
            (river_cells_y - gauge_y_idx) ** 2 + (river_cells_x - gauge_x_idx) ** 2
        ).argmin()

        # the cell cordinates
        y, x = grid.coordinatesOf(river_cells_y[nn], river_cells_x[nn])
        return Gauge(id=gauge.id, y=y, x=x, size=size)
