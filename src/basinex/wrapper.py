# -*- coding: utf-8 -*-
# from collections import namedtuple


class NcFile(object):
    def __init__(self, fname, ydim, xdim, outpath=None, y_shift=1, x_shift=0):
        self.fname = fname
        self.ydim = ydim
        self.xdim = xdim
        self.outpath = outpath
        self.y_shift = y_shift
        self.x_shift = x_shift


class GridFile(object):
    def __init__(self, fname, outpath=None):
        self.fname = fname
        self.outpath = outpath
