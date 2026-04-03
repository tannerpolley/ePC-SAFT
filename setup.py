import os

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup

try:
    from includeigen import get_include as get_eigen_include
except ImportError as exc:
    raise RuntimeError(
        "pcsaft now expects Eigen headers from the installed includeigen package. "
        "Install it first with `python -m pip install includeigen`."
    ) from exc


PACKAGE_ROOT = "src/pcsaft"

extra_compile_args = []
if os.name == "nt":
    extra_compile_args.append("/wd4551")

ext_modules = [
    Extension(
        "pcsaft.pcsaft",
        sources=[
            f"{PACKAGE_ROOT}/pcsaft.pyx",
            f"{PACKAGE_ROOT}/pcsaft_electrolyte.cpp",
        ],
        language="c++",
        include_dirs=[
            np.get_include(),
            PACKAGE_ROOT,
            get_eigen_include(),
        ],
        extra_compile_args=extra_compile_args,
    )
]

setup(
    ext_modules=cythonize(
        ext_modules,
        language_level="3",
    )
)
