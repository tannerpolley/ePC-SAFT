"""PEP 517 backend wrapper for sandbox-safe Windows package builds."""

from __future__ import annotations

import errno
import os
import shutil
import sys
import tempfile
from pathlib import Path

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
    if not os.environ.get("CMAKE_GENERATOR") and not shutil.which("cl"):
        if shutil.which("ninja"):
            os.environ.setdefault("CMAKE_GENERATOR", "Ninja")
        elif shutil.which("mingw32-make"):
            os.environ.setdefault("CMAKE_GENERATOR", "MinGW Makefiles")
            os.environ.setdefault("CMAKE_MAKE_PROGRAM", shutil.which("mingw32-make") or "")


def _has_build_dir(config_settings) -> bool:
    if not config_settings:
        return False
    return any(str(key).replace("_", "-") == "build-dir" for key in config_settings)


def _source_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _is_under(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def _external_temp_root() -> Path | None:
    source_root = _source_root()
    candidates: list[Path] = []
    local_app_data = os.environ.get("LOCALAPPDATA")
    if local_app_data:
        candidates.append(Path(local_app_data) / "Temp")
    if os.name != "nt":
        candidates.append(Path("/tmp"))
    candidates.append(Path(tempfile.gettempdir()))

    for candidate in candidates:
        try:
            resolved = candidate.expanduser().resolve()
            if _is_under(resolved, source_root):
                continue
            resolved.mkdir(parents=True, exist_ok=True)
        except OSError:
            continue
        return resolved
    return None


def _isolated_build_config(config_settings=None):
    config = dict(config_settings or {})
    if _has_build_dir(config):
        return config
    persistent = os.environ.get("EPCSAFT_PEP517_BUILD_DIR")
    if persistent:
        build_dir = Path(persistent).expanduser().resolve()
        build_dir.mkdir(parents=True, exist_ok=True)
    else:
        build_dir = Path(tempfile.mkdtemp(prefix="epcsaft-pep517-build-", dir=_external_temp_root())).resolve()
    config["build-dir"] = str(build_dir)
    return config


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    return _scikit_build.build_wheel(
        wheel_directory,
        config_settings=_isolated_build_config(config_settings),
        metadata_directory=metadata_directory,
    )


build_sdist = _scikit_build.build_sdist
get_requires_for_build_sdist = _scikit_build.get_requires_for_build_sdist
get_requires_for_build_wheel = _scikit_build.get_requires_for_build_wheel
prepare_metadata_for_build_wheel = _scikit_build.prepare_metadata_for_build_wheel
