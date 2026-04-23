import sys
from pathlib import Path

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_config

try:
    from includeigen import get_include as get_eigen_include
except ImportError as exc:
    raise RuntimeError(
        "epcsaft now expects Eigen headers from the installed includeigen package. "
        "Install it first with `python -m pip install includeigen`."
    ) from exc


ext_modules = [
    Extension(
        "epcsaft.epcsaft",
        sources=build_config.extension_sources(),
        language="c++",
        include_dirs=build_config.include_dirs(np.get_include(), get_eigen_include()),
        extra_compile_args=build_config.extra_compile_args(),
    )
]

setup(
    ext_modules=cythonize(
        ext_modules,
        language_level="3",
    )
)

