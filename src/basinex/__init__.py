# -*- coding: utf-8 -*-
"""
Purpose
=======

basinex is a basin extractor.
"""

try:
    from ._version import __version__
except ModuleNotFoundError:  # pragma: nocover
    # package is not installed
    __version__ = "0.0.0.dev0"

__all__ = ["__version__"]
