#cython: language_level=3
# distutils: language = c++
# -*- coding: utf-8 -*-

cimport cython

import numpy as np

cimport numpy as np
from libcpp.pair cimport pair
from libcpp.stack cimport stack
from libcpp.vector cimport vector

ctypedef pair[long, long] index

cdef vector[int] Y_OFFSET = [-1, -1, -1,  0,  0,   1,  1,  1]
cdef vector[int] X_OFFSET = [-1,  0,  1, -1,  1,  -1,  0,  1]
cdef vector[int] UPSTREAM_FDIRS = [ 2,  4,  8,  1, 16, 128, 64, 32]


cpdef char[:,:] extract(int[:,:] fdir, long gauge_y, long gauge_x):
    cdef char[:,:] mask = np.zeros_like(fdir, dtype=np.int8)
    cdef stack[index] stck
    cdef index idx

    stck.push(index(gauge_y, gauge_x))
    while (not stck.empty()):
        idx = stck.top()
        stck.pop()
        mask[idx.first, idx.second] = 1
        upstreamCells(fdir, stck, idx)

    return mask

@cython.boundscheck(False)
cdef void upstreamCells(int[:, :]& fdir, stack[index]& stack, index& idx):
    cdef int nrows = fdir.shape[0]
    cdef int ncols = fdir.shape[1]
    cdef long ynn, xnn
    cdef int i
    for i in range(8):
        ynn = idx.first + Y_OFFSET[i]
        xnn = idx.second + X_OFFSET[i]
        if (ynn >= 0) and (ynn < nrows) and (xnn >= 0) and (xnn < ncols):
            if fdir[ynn, xnn] == UPSTREAM_FDIRS[i]:
                stack.push(index(ynn, xnn))



