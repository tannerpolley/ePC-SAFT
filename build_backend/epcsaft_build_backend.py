"""PEP 517 backend wrapper for sandbox-safe Windows package builds."""

from __future__ import annotations

import errno
import os
import shutil
import sys
import tempfile

from scikit_build_core import build as _scikit_build


def _sandbox_safe_mkdtemp(suffix=None, prefix=None, dir=None):
    prefix, suffix, dir, output_type = tempfile._sanitize_params(prefix, suffix, dir)

    names = tempfile._get_candidate_names()
    if output_type is bytes:
        names = map(os.fsencode, names)

    for _ in range(tempfile.TMP_MAX):
        name = next(names)
        path = os.path.join(dir, prefix + name + suffix)
        sys.audit("tempfile.mkdtemp", path)
        try:
            os.mkdir(path, 0o777)
        except FileExistsError:
            continue
        except PermissionError:
            if os.name == "nt" and os.path.isdir(dir) and os.access(dir, os.W_OK):
                continue
            raise
        return os.path.abspath(path)

    raise FileExistsError(errno.EEXIST, "No usable temporary directory name found")


if os.name == "nt":
    tempfile.mkdtemp = _sandbox_safe_mkdtemp
    if shutil.which("mingw32-make") and not shutil.which("cl"):
        os.environ.setdefault("CMAKE_GENERATOR", "MinGW Makefiles")
        os.environ.setdefault("CMAKE_MAKE_PROGRAM", shutil.which("mingw32-make") or "")


build_sdist = _scikit_build.build_sdist
build_wheel = _scikit_build.build_wheel
build_editable = _scikit_build.build_editable
get_requires_for_build_sdist = _scikit_build.get_requires_for_build_sdist
get_requires_for_build_wheel = _scikit_build.get_requires_for_build_wheel
get_requires_for_build_editable = _scikit_build.get_requires_for_build_editable
prepare_metadata_for_build_wheel = _scikit_build.prepare_metadata_for_build_wheel
prepare_metadata_for_build_editable = _scikit_build.prepare_metadata_for_build_editable
