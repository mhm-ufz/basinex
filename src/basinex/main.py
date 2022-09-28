# -*- coding: utf-8 -*-

import logging
import os
import warnings
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import yaml

from . import __version__
from . import geoarray as ga
from .extractor import extract
from .gauges import matchFlowacc, readGauges
from .netcdf import NcDimDataset
from .netcdf4 import NcDataset
from .wrapper import GridFile, NcFile


def cli():
    parser = initArgparser()
    args = parser.parse_args()

    logging.basicConfig(
        format="%(message)s", level=logging.DEBUG if args.verbose else logging.INFO
    )

    os.chdir(Path(args.cwd))
    logging.debug(f"in directory: {os.getcwd()}")

    config = readConfig(args.input)

    # command line option -n
    line = args.line
    gauges = readGauges(
        config["gauges"], lat_fix=config.get("latitude-size-correction", False)
    )
    if line:
        if line > len(gauges) - 1:
            raise ValueError("Given line number exceeds table row count")
        gauges = gauges[line : line + 1]

    main(config, gauges)


def readConfig(fname):
    with open(fname, "r") as f:
        try:
            out = yaml.safe_load(f)
        except IOError:
            raise
        except Exception:
            raise RuntimeError("Syntactic error in {:}".format(fname))
    return out


def gaugeBasinMask(flowdir, gauge):
    gauge_idx = flowdir.indexOf(gauge.y, gauge.x)

    mask = np.asarray(
        extract(np.array(flowdir, dtype=np.int32, copy=True), *gauge_idx),
        dtype=np.int32,
    )

    mask[mask == 0] = flowdir.fill_value
    out = ga.array(mask, **flowdir.header)
    return out.trim()


def gridBasinMask(gauge):

    with NcDataset(gauge.path) as ncbase:
        # Infer dimension/fill_value information from the file
        var = ncbase.variables[gauge.varname]
        if var.ndim > 2:
            raise RuntimeError("Expected mask should not have more than 2 dimensions")

        y, x = var.dimensions
        # NOTE: would be a good to have y_shift and x_shift
        #        as optional fields in the gauge lut.
        with NcDimDataset(gauge.path, ydim=y, xdim=x, y_shift=0.5, x_shift=0.5) as nc:
            origin = nc.origin
            out = ga.array(
                nc.variables[gauge.varname][:],
                yorigin=nc.bbox["ymin" if origin[0] == "l" else "ymax"],
                xorigin=nc.bbox["xmin" if origin[1] == "l" else "xmax"],
                origin=origin,
                fill_value=var.fill_value,
                cellsize=nc.cellsize,
            )
    return out.setMask(out <= 0)


def openNcFiles(flist):
    out = {}
    flist = flist or ()
    for fdict in flist:
        fitem = NcFile(**fdict)
        out[fitem] = NcDimDataset(
            fitem.fname,
            "r",
            fitem.ydim,
            fitem.xdim,
            y_shift=fitem.y_shift,
            x_shift=fitem.x_shift,
        )

    return out


def openGridFiles(flist):
    out = {}
    flist = flist or ()
    for fdict in flist:
        out[GridFile(**fdict)] = ga.fromfile(fdict["fname"])
    return out


def writeFiles(bpath, fdict):
    for fitem, fobj in fdict.items():
        path = os.path.join(bpath, fitem.outpath or "")
        if not os.path.isdir(path):
            os.makedirs(path)
        fname = os.path.join(path, os.path.split(fitem.fname)[-1])
        logging.debug("writing file: %s", fname)
        fobj.tofile(fname)


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


def gaugeGrid(grid_template, gauge, out_in=None):
    if out_in is None:
        out = ga.full_like(grid_template, grid_template.fill_value)
        out._fobj = None
    else:
        out = out_in
    idx = out.indexOf(gauge.y, gauge.x)
    if out_in is not None and not out.mask[idx]:
        warnings.warn(f"There is already a gauge at {idx}, with id '{out.data[idx]}', replacing it by '{gauge.id}'.")
    out.data[idx] = gauge.id
    out.mask[idx] = False
    return out


def sameExtent(fobjs):
    bbox = commonBbox(fobjs)
    for fobj in fobjs:
        if fobj.bbox != bbox:
            return False
    return True


def writeReport(bpath, updated_gauge, gauge, error):
    Path(bpath).mkdir(exist_ok=True, parents=True)
    with open(os.path.join(bpath, "report.out"), "w") as f:
        f.write("error_catchment_size (%) : {:}\n".format(error))
        if updated_gauge is not None:
            f.write("new_catchment_size       : {:}\n".format(updated_gauge.size))
            f.write("new_y                    : {:}\n".format(updated_gauge.y))
            f.write("new_x                    : {:}\n".format(updated_gauge.x))
        f.write("input_catchment_size     : {:}\n".format(gauge.size))
        f.write("input_y                  : {:}\n".format(gauge.y))
        f.write("input_x                  : {:}\n".format(gauge.x))


def maskData(data, mask):
    if all(np.isclose(x, y) for x, y in zip(mask.cellsize, data.cellsize)):
        return data.setMask(mask.mask)

    enlarged_mask = mask.enlarge(**data.bbox).astype(float)
    rescaled_mask = ga.rescale(
        enlarged_mask, abs(data.cellsize[0] / mask.cellsize[0]), func="average"
    ).copy()
    return data.setMask(rescaled_mask.mask)


def main(config, gauges):

    flowacc = None
    flowdir = None
    filedict_main = {}
    gaugefile_main = None
    updated_gauge = None

    for gauge in gauges:
        logging.info("processing gauge: %s", gauge.id)

        filedict = {}

        if not gauge.path:

            # create mask if not given
            if flowacc is None:
                logging.debug("reading flow accumulation")
                flowacc = ga.fromfile(config["flowacc"])

            if flowdir is None:
                logging.debug("reading flow direction")
                flowdir = ga.fromfile(config["flowdir"]).astype(np.int32)

            if gauge.size:
                logging.debug("moving gauge to streamflow")
                updated_gauge, error = matchFlowacc(gauge, flowacc, **config["matching"])

            if not updated_gauge:
                warnings.warn("Failed to match the gauge to the flow accumulation grid")
                continue

            logging.debug("generating basin mask")
            mask = gaugeBasinMask(flowdir, updated_gauge)

            # write gauge grid if desired
            if "gauge" in config:
                logging.debug("writing gauge file")
                fname = config["gauge"].get("fname", "idgauges.asc")
                if len(gauges) > 1:
                    if gaugefile_main is None:
                        fitem_main = GridFile(
                            fname=fname,
                            outpath=config["gauge"].get("outpath"),
                        )
                    gaugefile_main = gaugeGrid(flowacc, updated_gauge, gaugefile_main)
                    filedict_main[fitem_main] = gaugefile_main
                fitem = GridFile(
                    fname=fname,
                    outpath=config["gauge"].get("outpath"),
                )
                gaugefile = gaugeGrid(flowacc, updated_gauge).shrink(**mask.bbox)
                filedict[fitem] = maskData(gaugefile, mask)

        else:
            logging.debug("reading gauge file")
            mask = gridBasinMask(gauge)

        for fdict in config.get("gridfiles", []):
            logging.debug("processing: %s", fdict["fname"])
            griddata = ga.fromfile(fdict["fname"]).shrink(**mask.bbox)
            filedict[GridFile(**fdict)] = maskData(griddata, mask)

        for fdict in config.get("ncfiles", []):
            logging.debug("processing: %s", fdict["fname"])
            fitem = NcFile(**fdict)
            ncdata = NcDimDataset(
                fitem.fname,
                "r",
                fitem.ydim,
                fitem.xdim,
                y_shift=fitem.y_shift,
                x_shift=fitem.x_shift,
            )
            filedict[NcFile(**fdict)] = maskData(ncdata.shrink(**mask.bbox), mask)

        # write mask grid if desired
        if "mask" in config:
            logging.debug("writing mask file")
            fitem = GridFile(
                fname=config["mask"].get("fname", "mask.asc"),
                outpath=config["mask"].get("outpath"),
            )
            filedict[fitem] = mask

        if filedict:
            logging.debug("finding common extent")
            bbox = commonBbox(tuple(filedict.values()))

            logging.debug("enlarging data to common extent")
            filedict = enlargeFiles(filedict, bbox)

            if not sameExtent(tuple(filedict.values())):
                raise RuntimeError("incompatible cellsizes")

        bpath = os.path.join(config["outpath"], gauge.id)
        writeFiles(bpath, filedict)
        logging.debug("writing report")
        writeReport(bpath, updated_gauge, gauge, error)

    writeFiles(config["outpath"], filedict_main)


def initArgparser():

    parser = ArgumentParser(description="mHM basin extractor")

    parser.add_argument(
        "-n",
        "--line",
        type=int,
        help=(
            "the gauge to extract, given as its (0-based) "
            "line number in the look up table"
        ),
    )

    parser.add_argument(
        "-i",
        "--input",
        default="input.yml",
        help="the input yaml file to read (default: 'input.yml')",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="give some status output",
    )

    parser.add_argument("-c", "--cwd", default=".", help="the working directory")

    parser.add_argument("--version", action="version", version=__version__)

    return parser


if __name__ == "__main__":
    cli()
