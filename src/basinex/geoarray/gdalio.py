# -*- coding: utf-8 -*-

import os
import warnings

import numpy as np

try:
    from osgeo import gdal, osr
except ImportError:
    import gdal, osr

from .gdalspatial import _Projection
from .geotrans import _Geotrans

gdal.UseExceptions()
gdal.PushErrorHandler("CPLQuietErrorHandler")

# should be extended, for available options see:
# http://www.gdal.org/formats_list.html
_DRIVER_DICT = {
    ".tif": "GTiff",
    ".asc": "AAIGrid",
    ".img": "HFA",
    ".sdat": "SAGA",
    ".png": "PNG",  # not working properly
}

# type mapping:
#     - there is no boolean data type in GDAL
_TYPEMAP = {
    "uint8": 1,
    "int8": 1,
    "uint16": 2,
    "int16": 3,
    "uint32": 4,
    "int32": 5,
    # "int64"      : 5, # there is no int64 data type in GDAL, map to int32 and issue a warning
    "float32": 6,
    "float64": 7,
    "complex64": 10,
    "complex128": 11,
    1: "int8",
    2: "uint16",
    3: "int16",
    4: "uint32",
    5: "int32",
    6: "float32",
    7: "float64",
    10: "complex64",
    11: "complex128",
}

_COLOR_DICT = {
    1: "L",
    2: "P",
    3: "R",
    4: "G",
    5: "B",
    6: "A",
    7: "H",
    8: "S",
    9: "V",
    10: "C",
    11: "M",
    12: "Y",
    13: "K",
    14: "Y",
    15: "Cb",
    16: "Cr",
}

_COLOR_MODE_LIST = ("L", "P", "RGB", "RGBA", "CMYK", "HSV", "YCbCr")

_FILE_MODE_DICT = {"r": gdal.GA_ReadOnly, "v": gdal.GA_ReadOnly, "a": gdal.GA_Update}


def _fromFile(fname, mode="r"):
    """
    Parameters
    ----------
    fname : str  # file name

    Returns
    -------
    GeoArray

    Purpose
    -------
    Create GeoArray from file

    """

    if mode not in _FILE_MODE_DICT:
        raise TypeError(
            "Supported file modes are: {:}".format(", ".join(_FILE_MODE_DICT.keys()))
        )

    fobj = gdal.OpenShared(fname, _FILE_MODE_DICT[mode])
    if fobj:
        return _fromDataset(fobj, mode)
    raise IOError("Could not open file: {:}".format(fname))


def _getColorMode(fobj):
    tmp = []
    for i in range(fobj.RasterCount):
        color = fobj.GetRasterBand(i + 1).GetColorInterpretation()
        tmp.append(_COLOR_DICT.get(color, "L"))
    return "".join(sorted(set(tmp), key=tmp.index))


def _fromDataset(fobj, mode="r"):

    from .core import GeoArray

    def _parseGeotrans(geotrans):
        return {
            "yorigin": geotrans[3],
            "xorigin": geotrans[0],
            "ycellsize": geotrans[5],
            "xcellsize": geotrans[1],
            "origin": "ul",  # is that always true?
            "yparam": geotrans[4],
            "xparam": geotrans[2],
        }

    fill_values = tuple(
        fobj.GetRasterBand(i + 1).GetNoDataValue() for i in range(fobj.RasterCount)
    )

    if len(set(fill_values)) > 1:
        warnings.warn(
            "More then on fill value found. Only {:} will be used".format(
                fill_values[0]
            ),
            RuntimeWarning,
        )

    data = fobj.GetVirtualMemArray() if mode == "v" else fobj.ReadAsArray()
    # NOTE: not to robust...
    geotrans = _Geotrans(shape=data.shape, **_parseGeotrans(fobj.GetGeoTransform()))

    return GeoArray(
        data=data,
        fill_value=fill_values[0],
        proj=_Projection(fobj.GetProjection()),
        mode=mode,
        color_mode=_getColorMode(fobj),
        fobj=fobj,
        geotrans=geotrans,
    )


def _getDataset(grid, mem=False):

    # Returns an gdal memory dataset created from the given grid

    if grid._fobj and not mem:
        return grid._fobj

    driver = gdal.GetDriverByName("MEM")

    try:
        out = driver.Create(
            "", grid.ncols, grid.nrows, grid.nbands, _TYPEMAP[str(grid.dtype)]
        )
    except KeyError:
        raise RuntimeError("Datatype {:} not supported by GDAL".format(grid.dtype))

    if isinstance(grid.geotrans, _Geotrans):
        out.SetGeoTransform(grid.toGdal())

    if grid.proj:
        out.SetProjection(grid.proj.toWkt())

    for n in range(grid.nbands):
        band = out.GetRasterBand(n + 1)
        if grid.fill_value is not None:
            band.SetNoDataValue(float(grid.fill_value))
        data = grid[n] if grid.ndim > 2 else grid
        band.WriteArray(data)

    return out


def _toFile(geoarray, fname):
    """
    Arguments
    ---------
    fname : str  # file name

    Returns
    -------
    None

    Purpose
    -------
    Write GeoArray to file. The output dataset type is derived from
    the file name extension. See _DRIVER_DICT for implemented formats.
    """

    def _fnameExtension(fname):
        return os.path.splitext(fname)[-1].lower()

    def _getDriver(fext):
        """
        Guess driver from file name extension
        """
        if fext in _DRIVER_DICT:
            driver = gdal.GetDriverByName(_DRIVER_DICT[fext])
            metadata = driver.GetMetadata_Dict()
            if "YES" == metadata.get("DCAP_CREATE", metadata.get("DCAP_CREATECOPY")):
                return driver
            raise IOError("Datatype cannot be written")
        raise IOError("No driver found for filename extension '{:}'".format(fext))

    def _getDatatype(driver):
        tnames = tuple(driver.GetMetadata_Dict()["DMD_CREATIONDATATYPES"].split(" "))
        types = tuple(gdal.GetDataTypeByName(t) for t in tnames)
        tdict = tuple((gdal.GetDataTypeSize(t), t) for t in types)
        otype = max(tdict, key=lambda x: x[0])[-1]
        return np.dtype(_TYPEMAP[otype])

    dataset = _getDataset(geoarray)
    driver = _getDriver(_fnameExtension(fname))
    driver.CreateCopy(fname, dataset, 0)


def _writeData(grid):

    fobj = grid.fobj
    data = grid.data
    if data.ndim == 2:
        data = data[None, ...]

    for n in range(fobj.RasterCount):
        fobj.GetRasterBand(n + 1).WriteArray(data[n])
