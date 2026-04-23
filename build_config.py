"""Shared build metadata for the ePC-SAFT Cython/C++ extension."""

from __future__ import annotations

import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
PACKAGE_ROOT = REPO_ROOT / "src" / "epcsaft"
NATIVE_ROOT = PACKAGE_ROOT / "native"
BUILD_ROOT = REPO_ROOT / "build"
BUILD_STAMP = BUILD_ROOT / "epcsaft-editable-build.json"

PACKAGE_ROOT_REL = "src/epcsaft"
NATIVE_ROOT_REL = f"{PACKAGE_ROOT_REL}/native"
CONTRIB_ROOT_REL = f"{NATIVE_ROOT_REL}/contributions"

CYTHON_SOURCE = f"{PACKAGE_ROOT_REL}/epcsaft.pyx"
PYX_SUPPORT = (
    f"{PACKAGE_ROOT_REL}/epcsaft.pxd",
)
NATIVE_HEADERS = (
    f"{NATIVE_ROOT_REL}/epcsaft_electrolyte.h",
    f"{NATIVE_ROOT_REL}/epcsaft_core_internal.h",
    f"{NATIVE_ROOT_REL}/epcsaft_autodiff_internal.h",
    f"{CONTRIB_ROOT_REL}/epcsaft_contrib_internal.h",
)
NATIVE_SOURCES = (
    f"{NATIVE_ROOT_REL}/epcsaft_parameter_setup.cpp",
    f"{NATIVE_ROOT_REL}/epcsaft_density.cpp",
    f"{CONTRIB_ROOT_REL}/epcsaft_contrib_hc.cpp",
    f"{CONTRIB_ROOT_REL}/epcsaft_contrib_disp.cpp",
    f"{CONTRIB_ROOT_REL}/epcsaft_contrib_assoc.cpp",
    f"{CONTRIB_ROOT_REL}/epcsaft_contrib_ion.cpp",
    f"{CONTRIB_ROOT_REL}/epcsaft_contrib_born.cpp",
    f"{NATIVE_ROOT_REL}/epcsaft_ares.cpp",
    f"{NATIVE_ROOT_REL}/epcsaft_thermo.cpp",
    f"{NATIVE_ROOT_REL}/epcsaft_Z.cpp",
    f"{NATIVE_ROOT_REL}/epcsaft_mu.cpp",
    f"{NATIVE_ROOT_REL}/epcsaft_fugcoef.cpp",
    f"{NATIVE_ROOT_REL}/epcsaft_activity.cpp",
    f"{NATIVE_ROOT_REL}/epcsaft_state.cpp",
    f"{NATIVE_ROOT_REL}/epcsaft_regression.cpp",
)
EXTENSION_SOURCES = (CYTHON_SOURCE, *NATIVE_SOURCES)
REBUILD_INPUTS = (
    "pyproject.toml",
    "setup.py",
    "MANIFEST.in",
    "build_config.py",
    CYTHON_SOURCE,
    *PYX_SUPPORT,
    *NATIVE_HEADERS,
    *NATIVE_SOURCES,
)


def repo_path(relative_path: str) -> Path:
    return REPO_ROOT / relative_path


def rebuild_input_paths() -> list[Path]:
    return [repo_path(path) for path in REBUILD_INPUTS]


def extension_sources() -> list[str]:
    return list(EXTENSION_SOURCES)


def include_dirs(numpy_include: str, eigen_include: str) -> list[str]:
    return [
        numpy_include,
        PACKAGE_ROOT_REL,
        NATIVE_ROOT_REL,
        CONTRIB_ROOT_REL,
        eigen_include,
    ]


def extra_compile_args() -> list[str]:
    if os.name == "nt":
        return ["/std:c++17", "/wd4551"]
    return ["-std=c++17"]


def generated_cython_cpp() -> Path:
    return PACKAGE_ROOT / "epcsaft.cpp"
