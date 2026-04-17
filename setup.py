import os
import sys
from pathlib import Path

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


def _ipopt_build_config():
    prefix = Path(os.environ.get("CONDA_PREFIX", sys.prefix)).resolve()
    if os.name == "nt":
        include_root = prefix / "Library" / "include"
        coin_include = include_root / "coin"
        lib_root = prefix / "Library" / "lib"
        header = coin_include / "IpIpoptApplication.hpp"
        lib = lib_root / "ipopt.lib"
    else:
        include_root = prefix / "include"
        coin_include = include_root / "coin"
        lib_root = prefix / "lib"
        header = coin_include / "IpIpoptApplication.hpp"
        lib = lib_root / "libipopt.so"

    if not header.exists() or not lib.exists():
        raise RuntimeError(
            "ePC-SAFT now requires IPOPT headers and libraries in the active environment. "
            f"Expected {header} and {lib}. Install IPOPT in the ePC-SAFT conda env first."
        )

    return {
        "include_dirs": [str(include_root), str(coin_include)],
        "library_dirs": [str(lib_root)],
        "libraries": ["ipopt"],
    }


ipopt_build = _ipopt_build_config()

extra_compile_args = []
if os.name == "nt":
    extra_compile_args.append("/wd4551")

ext_modules = [
    Extension(
        "epcsaft.epcsaft",
        sources=[
            f"{PACKAGE_ROOT}/epcsaft.pyx",
            f"{NATIVE_ROOT}/epcsaft_parameter_setup.cpp",
            f"{NATIVE_ROOT}/epcsaft_density.cpp",
            f"{NATIVE_ROOT}/contributions/epcsaft_contrib_hc.cpp",
            f"{NATIVE_ROOT}/contributions/epcsaft_contrib_disp.cpp",
            f"{NATIVE_ROOT}/contributions/epcsaft_contrib_assoc.cpp",
            f"{NATIVE_ROOT}/contributions/epcsaft_contrib_ion.cpp",
            f"{NATIVE_ROOT}/contributions/epcsaft_contrib_born.cpp",
            f"{NATIVE_ROOT}/epcsaft_ares.cpp",
            f"{NATIVE_ROOT}/epcsaft_thermo.cpp",
            f"{NATIVE_ROOT}/epcsaft_Z.cpp",
            f"{NATIVE_ROOT}/epcsaft_mu.cpp",
            f"{NATIVE_ROOT}/epcsaft_fugcoef.cpp",
            f"{NATIVE_ROOT}/epcsaft_activity.cpp",
            f"{NATIVE_ROOT}/epcsaft_state.cpp",
            f"{NATIVE_ROOT}/epcsaft_regression.cpp",
        ],
        language="c++",
        include_dirs=[
            np.get_include(),
            PACKAGE_ROOT,
            NATIVE_ROOT,
            f"{NATIVE_ROOT}/contributions",
            get_eigen_include(),
            *ipopt_build["include_dirs"],
        ],
        library_dirs=ipopt_build["library_dirs"],
        libraries=ipopt_build["libraries"],
        extra_compile_args=extra_compile_args,
    )
]

setup(
    ext_modules=cythonize(
        ext_modules,
        language_level="3",
    )
)

