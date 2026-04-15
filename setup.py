import os

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup

try:
    from includeigen import get_include as get_eigen_include
except ImportError as exc:
    raise RuntimeError(
        "epcsaft now expects Eigen headers from the installed includeigen package. "
        "Install it first with `python -m pip install includeigen`."
    ) from exc


PACKAGE_ROOT = "src/epcsaft"
NATIVE_ROOT = f"{PACKAGE_ROOT}/native"

extra_compile_args = []
if os.name == "nt":
    extra_compile_args.append("/wd4551")

ext_modules = [
    Extension(
        "epcsaft.epcsaft",
        sources=[
            f"{PACKAGE_ROOT}/epcsaft.pyx",
            f"{NATIVE_ROOT}/epcsaft_support_common.cpp",
            f"{NATIVE_ROOT}/epcsaft_support_assoc.cpp",
            f"{NATIVE_ROOT}/epcsaft_support_electrolyte.cpp",
            f"{NATIVE_ROOT}/epcsaft_ares.cpp",
            f"{NATIVE_ROOT}/epcsaft_dadrho.cpp",
            f"{NATIVE_ROOT}/epcsaft_dadx.cpp",
            f"{NATIVE_ROOT}/epcsaft_dadt.cpp",
            f"{NATIVE_ROOT}/epcsaft_Z.cpp",
            f"{NATIVE_ROOT}/epcsaft_mu.cpp",
            f"{NATIVE_ROOT}/epcsaft_fugcoef.cpp",
            f"{NATIVE_ROOT}/epcsaft_activity.cpp",
            f"{NATIVE_ROOT}/epcsaft_density_solver.cpp",
            f"{NATIVE_ROOT}/epcsaft_state.cpp",
        ],
        language="c++",
        include_dirs=[
            np.get_include(),
            PACKAGE_ROOT,
            NATIVE_ROOT,
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

