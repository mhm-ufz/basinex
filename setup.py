# -*- coding: utf-8 -*-
"""Basin Extractor."""
import os
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

# cython extensions ###########################################################

CY_MODULES = []
CY_MODULES.append(
    Extension(
        "basinex.extractor",
        [os.path.join("src", "basinex", "extractor.pyx")],
        language="c++",
        include_dirs=[np.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    )
)
EXT_MODULES = cythonize(CY_MODULES)  # annotate=True

# embed signatures for sphinx
for ext_m in EXT_MODULES:
    ext_m.cython_directives = {"embedsignature": True}

# setup #######################################################################

setup(ext_modules=EXT_MODULES, include_dirs=[np.get_include()])
