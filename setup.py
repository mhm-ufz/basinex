#! /usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension

from Cython.Build import cythonize
import numpy as np

# usage: python setup.py build_ext --inplace

extension = Extension(
    "extractor",
    sources=["extractor/extractor.pyx"],
    language="c++",
    include_dirs=[np.get_include()])

setup(ext_modules=cythonize(extension))

