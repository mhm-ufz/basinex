#! /usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup
from Cython.Build import cythonize

# usage: python setup.py build_ext --inplace

setup(
    ext_modules = cythonize("extractor/extractor.pyx",
                            language="c++",))

