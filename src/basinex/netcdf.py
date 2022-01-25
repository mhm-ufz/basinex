# -*- coding: utf-8 -*-

import os

import numpy as np

from .netcdf4 import NcDataset


class NcDimDataset(NcDataset):
    def __init__(
        self,
        fname,
        mode="r",
        ydim="y",
        xdim="x",
        cellsize=None,
        y_shift=1,
        x_shift=0,
        *args,
        **kwargs,
    ):
        # shift_factor: the fraction of a cellsize the origin is shifted
        super(NcDimDataset, self).__init__(fname, mode, *args, **kwargs)
        self.__dict__["_y"] = ydim
        self.__dict__["_x"] = xdim
        self.__dict__["_cellsize"] = cellsize
        self.__dict__["y_shift"] = float(y_shift)
        self.__dict__["x_shift"] = float(x_shift)

    # def _delta(self):
    # ycs, xcs = np.abs(self.cellsize)
    # return ycs * self.y_shift, xcs * self.x_shift

    @property
    def bbox(self):
        x = self.variables[self._x][:]
        y = self.variables[self._y][:]
        ycs, xcs = [abs(v) for v in self.cellsize]

        return {
            "ymin": np.min(y) - (ycs * (1 - self.y_shift)),
            "ymax": np.max(y) + (ycs * self.y_shift),
            "xmin": np.min(x) - (xcs * (1 - self.x_shift)),
            "xmax": np.max(x) + (xcs * self.x_shift),
        }

    @property
    def origin(self):
        dy, dx = self.cellsize
        origin = ("l" if dy >= 0 else "u", "l" if dx >= 0 else "r")
        return "".join(origin)

    @property
    def cellsize(self):
        if self._cellsize is None:
            y = self.variables[self._y][:2]
            x = self.variables[self._x][:2]
            try:
                self.__dict__["_cellsize"] = (y[1] - y[0], x[1] - x[0])
            except IndexError:
                raise RuntimeError("Not enough coordinate values to infer cellsize")
        return self._cellsize

    @property
    def fill_values(self):
        return {k: v.fill_value for k, v in self.variables.items()}

    def shrink(self, ymin, ymax, xmin, xmax):
        def _removeCells(ymin, ymax, xmin, xmax):
            cellsize = [float(abs(v)) for v in self.cellsize]
            out = {
                "top": int(np.floor((self.bbox["ymax"] - ymax) / cellsize[0])),
                "left": int(np.floor((xmin - self.bbox["xmin"]) / cellsize[1])),
                "bottom": int(np.floor((ymin - self.bbox["ymin"]) / cellsize[0])),
                "right": int(np.floor((self.bbox["xmax"] - xmax) / cellsize[1])),
            }
            return out

        def _slices(cells):

            y = self.variables[self._y][:]
            ystart, yend = (
                (cells["bottom"], cells["top"])
                if y[0] <= y[-1]
                else (cells["top"], cells["bottom"])
            )
            yend = len(y) - yend

            x = self.variables[self._x][:]
            xstart, xend = (
                (cells["left"], cells["right"])
                if x[0] <= x[-1]
                else (cells["right"], cells["left"])
            )
            xend = len(x) - xend

            return {
                self._y: slice(ystart, yend),
                self._x: slice(xstart, xend),
            }

        cells = _removeCells(ymin, ymax, xmin, xmax)
        shrink_slices = _slices(cells)

        nc = NcDimDataset(
            None, "w", self._y, self._x, self.cellsize, self.y_shift, self.x_shift
        )
        nc.copyAttributes(self.attributes)
        nc.copyDimensions(self.dimensions, skip=(self._y, self._x))
        nc.createDimensions({k: (v.stop - v.start) for k, v in shrink_slices.items()})

        for name, var in self.variables.items():
            if any(d in var.dimensions for d in shrink_slices.keys()):
                slices = [shrink_slices.get(d, slice(None)) for d in var.dimensions]
                newvar = nc.copyVariable(var, data=False)
                newvar[:] = var[slices]
            else:
                nc.copyVariable(var, data=True)
        return nc

    def enlarge(self, ymin, ymax, xmin, xmax):
        def _padCells(ymin, ymax, xmin, xmax):
            cellsize = [float(abs(v)) for v in self.cellsize]

            out = {
                "top": int(np.ceil((ymax - self.bbox["ymax"]) / cellsize[0])),
                "left": int(np.ceil((self.bbox["xmin"] - xmin) / cellsize[1])),
                "bottom": int(np.ceil((self.bbox["ymin"] - ymin) / cellsize[0])),
                "right": int(np.ceil((xmax - self.bbox["xmax"]) / cellsize[1])),
            }
            return out

        def _slices(padding):

            y = self.variables[self._y][:]
            ystart = padding["bottom" if y[0] <= y[-1] else "top"]
            yend = ystart + len(y)

            x = self.variables[self._x][:]
            xstart = padding["left" if x[0] <= x[-1] else "right"]
            xend = xstart + len(x)

            return {self._y: slice(ystart, yend), self._x: slice(xstart, xend)}

        def _coordinates(padding):
            ycs, xcs = self.cellsize

            yvals = self.variables[self._y][:]
            xvals = self.variables[self._x][:]

            ynum = len(self.dimensions[self._y]) + padding["bottom"] + padding["top"]
            ymin = np.min(yvals) - padding["bottom"] * abs(ycs)
            ymax = np.max(yvals) + padding["bottom"] * abs(ycs)

            xnum = len(self.dimensions[self._x]) + padding["bottom"] + padding["top"]
            xmin = np.min(xvals) - padding["bottom"] * abs(xcs)
            xmax = np.max(xvals) + padding["bottom"] * abs(xcs)

            out = {
                self._y: sorted(
                    np.linspace(ymin, ymax, num=ynum, endpoint=True), reverse=ycs < 0
                ),
                self._x: sorted(
                    np.linspace(xmin, xmax, num=xnum, endpoint=True), reverse=xcs < 0
                ),
            }
            return out

        padding = _padCells(ymin, ymax, xmin, xmax)
        enlarge_slices = _slices(padding)
        coordinates = _coordinates(padding)

        nc = NcDimDataset(
            None, "w", self._y, self._x, self.cellsize, self.y_shift, self.x_shift
        )

        nc.copyAttributes(self.attributes)
        nc.copyDimensions(self.dimensions, skip=(self._y, self._x))
        nc.createDimensions({k: (len(v)) for k, v in coordinates.items()})

        for name, var in self.variables.items():
            if any(d in var.dimensions for d in enlarge_slices.keys()):
                slices = [enlarge_slices.get(d, slice(None)) for d in var.dimensions]
                newvar = nc.copyVariable(var, data=False)
                if name in (self._y, self._x):
                    newvar[:] = coordinates[name]
                else:
                    newvar[slices] = var[:]
            else:
                nc.copyVariable(var, data=True)

        return nc

    def setMask(self, mask):
        assert mask.ndim == 2, "expected a 2D mask array"

        y_idx, x_idx = np.nonzero(mask)

        mask_slices = {
            self._y: y_idx,
            self._x: x_idx,
        }

        nc = NcDimDataset(
            None, "w", self._y, self._x, self.cellsize, self.y_shift, self.x_shift
        )
        nc.copyAttributes(self.attributes)
        nc.copyDimensions(self.dimensions)

        for name, var in self.variables.items():
            newvar = nc.copyVariable(var, data=True)
            if self._y in var.dimensions and self._x in var.dimensions:
                slices = [mask_slices.get(d, slice(None)) for d in var.dimensions]
                # do it the clunky way because of a slicing bug in netCDF4
                data = newvar[:]
                data[tuple(slices)] = var.fill_value
                newvar[:] = data
        return nc

    def _header(self):
        fill_value = tuple(v for v in self.fill_values.values() if v) or (None,)
        return "\n".join(
            [
                "ncols\t{0}".format(len(self.variables[self._x])),
                "nrows\t{0}".format(len(self.variables[self._y])),
                "xllcorner\t{0}".format(self.bbox["xmin"]),
                "yllcorner\t{0}".format(self.bbox["ymin"]),
                "cellsize\t{0}".format(abs(self.cellsize[0])),
                "NODATA_value\t{0}\n".format(fill_value[0]),
            ]
        )

    def tofile(self, fname):
        with NcDataset(fname, "w") as nc:
            nc.copyDataset(self, vardata=True)
        # write a header file
        path = os.path.split(fname)[0]
        with open(os.path.join(path, "header.txt"), "w") as f:
            f.write(self._header())
