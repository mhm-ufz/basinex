#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
from argparse import ArgumentParser
import yaml
import numpy as np
from ufz.netcdf4 import NcDataset
from lib import geoarray as ga
from lib.netcdf import NcDimDataset
from lib.gauges import readGauges, matchFlowacc
from lib.wrapper import NcFile, GridFile
from extractor import extract
import warnings


INPUT = "input.yml"


def readConfig(fname):
    with open(fname, "r") as f:
        try:
            out = yaml.load(f)
        except IOError:
            raise
        except Exception:
            raise RuntimeError("Syntactic error in {:}".format(fname))
    return out

        
def gaugeBasinMask(flowdir, gauge):
    gauge_idx = flowdir.indexOf(gauge.y, gauge.x)

    mask = np.asarray(
        extract(
            np.array(flowdir, dtype=np.int32, copy=True),
            *gauge_idx),
        dtype=np.int32)

    mask[mask == 0] = flowdir.fill_value
    out = ga.array(mask, **flowdir.header)
    return out.trim()


def gridBasinMask(gauge):

    with NcDataset(gauge.path) as ncbase:
        # Infer dimension/fill_value information from the file
        var = ncbase.variables[gauge.varname]
        if var.ndim > 2:
            raise RuntimeError(
                "Expected mask should not have more than 2 dimensions")

        y, x = var.dimensions
        with NcDimDataset(gauge.path, ydim=y, xdim=x) as nc:
            origin = nc.origin
            out = ga.array(
                nc.variables[gauge.varname][:],
                yorigin=nc.bbox["ymin" if origin[0] == "l" else "ymax"],
                xorigin=nc.bbox["xmin" if origin[1] == "l" else "xmax"],
                origin=origin,
                fill_value=var.fill_value,
                cellsize=nc.cellsize
            )
    return out.setMask(out <= 0)


def openNcFiles(flist):
    out = {}
    for fdict in flist:
        fitem = NcFile(**fdict)
        out[fitem] = NcDimDataset(
            fitem.fname, "r", fitem.ydim, fitem.xdim,
            y_shift=fitem.y_shift, x_shift=fitem.x_shift)

    return out


def openGridFiles(flist):
    out = {}
    for fdict in flist:
        out[GridFile(**fdict)] = ga.fromfile(fdict["fname"])
    return out


def writeFiles(bpath, fdict):
    for fitem, fobj in fdict.items():
        path = os.path.join(bpath, fitem.outpath or "")
        if not os.path.isdir(path):
            os.makedirs(path)
        fobj.tofile(os.path.join(path, os.path.split(fitem.fname)[-1]))


def enlargeFiles(fdict, bbox):
    out = {}
    for fitem, fobj in fdict.items():
        out[fitem] = fobj.enlarge(**bbox)

    return out


def commonBbox(fobjs):
    bbox = fobjs[0].bbox
    for fobj in fobjs[1:]:
        bbox["ymin"] = min(bbox["ymin"], fobj.bbox["ymin"])
        bbox["ymax"] = max(bbox["ymax"], fobj.bbox["ymax"])
        bbox["xmin"] = min(bbox["xmin"], fobj.bbox["xmin"])
        bbox["xmax"] = max(bbox["xmax"], fobj.bbox["xmax"])
    return bbox


def gaugeGrid(grid_template, gauge):
    out = ga.full_like(grid_template, grid_template.fill_value)
    idx = out.indexOf(gauge.y, gauge.x)
    out.data[idx] = gauge.id
    out.mask[idx] = False
    return out


def sameExtend(fobjs):
    bbox = commonBbox(fobjs)
    for fobj in fobjs:
        if fobj.bbox != bbox:
            return False
    return True


def writeReport(bpath, mask, scaling_factor):
    size = ((np.sum(~mask.mask) * np.prod(np.abs(mask.cellsize)))
            * scaling_factor**2)
    with open(os.path.join(bpath, "report.out"), "w") as f:
        f.write("calculated_catchment_size: {:}\n".format(size))


def maskData(data, mask):
    if all(x == y for x, y in zip(mask.cellsize, data.cellsize)):
        return data.setMask(mask.mask)

    enlarged_mask = mask.enlarge(**data.bbox).astype(float)
    rescaled_mask = ga.rescale(
        enlarged_mask, abs(data.cellsize[0]/mask.cellsize[0]),
        func='average').copy()
    return data.setMask(rescaled_mask.mask)


def main(config, gauges):

    for gauge in gauges:

        filedict={}

        if not gauge.path:

            # create mask if not given
            flowacc = ga.fromfile(config["flowacc"]).astype(np.int32)
            flowdir = ga.fromfile(config["flowdir"]).astype(np.int32)

            if gauge.size:
                # extract a single cell
                gauge = matchFlowacc(gauge, flowacc, **config["matching"])

            if not gauge:
                warnings.warn(
                    "Failed to match gauge '{:}' to the flowaccumulation grid"
                    .format(gauge.id))
                continue

            logging.debug("generating basin mask")
            mask = gaugeBasinMask(flowdir, gauge)

            # write gauge grid if desired
            if "gauge" in config:
                fitem = GridFile(
                    fname=config["gauge"].get("fname", "idgauges.asc"),
                    outpath=config["gauge"].get("outpath"))
                gaugefile = gaugeGrid(flowacc, gauge).shrink(**mask.bbox)
                filedict[fitem] = maskData(gaugefile, mask)

        else:
            mask = gridBasinMask(gauge)

        for fdict in config["gridfiles"]:
            logging.debug("processing: %s", fdict["fname"])
            griddata = ga.fromfile(fdict["fname"]).shrink(**mask.bbox)
            filedict[GridFile(**fdict)] = maskData(griddata, mask)

        for fdict in config["ncfiles"]:
            logging.debug("processing: %s", fdict["fname"])
            fitem = NcFile(**fdict)
            ncdata = NcDimDataset(
                fitem.fname, "r", fitem.ydim, fitem.xdim,
                y_shift=fitem.y_shift, x_shift=fitem.x_shift)
            filedict[NcFile(**fdict)] = maskData(ncdata.shrink(**mask.bbox), mask)

        # write mask grid if desired
        if "mask" in config:
            fitem = GridFile(
                fname=config["mask"].get("fname", "mask.asc"),
                outpath=config["mask"].get("outpath"))
            filedict[fitem] = mask

        bbox = commonBbox(filedict.values())

        filedict = enlargeFiles(filedict, bbox)

        if not sameExtend(filedict.values()):
            raise RuntimeError("incompatible cellsizes")

        bpath = os.path.join(config["outpath"], gauge.id)
        writeFiles(bpath, filedict)
        writeReport(bpath, mask, config["matching"]["scaling_factor"])
    

if __name__ == "__main__":

    parser = ArgumentParser(description="Extract basins")

    parser.add_argument(
        "-n", required=False, dest="line", type=int,
        help=("the gauge to extract, given as its (0-based) "
              "line number in the look up table"))

    args = parser.parse_args()
    line = args.line

    config = readConfig(INPUT)

    # command line option -n
    gauges = readGauges(config["gauges"])
    if line:
        if line > len(gauges)-1:
            raise ValueError("Given line number exceeds table row count")
        gauges = gauges[line:line+1]

    main(config, gauges)
