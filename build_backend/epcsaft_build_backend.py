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


def _truthy_env(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


def _config_has(config: dict, key: str) -> bool:
    return any(str(existing).replace("_", "-") == key.replace("_", "-") for existing in config)


def _set_config_default(config: dict, key: str, value: str) -> None:
    if not _config_has(config, key):
        config[key] = value


def _config_value(config: dict, key: str) -> str | None:
    normalized_key = key.replace("_", "-")
    for existing, value in config.items():
        if str(existing).replace("_", "-") == normalized_key:
            return str(value)
    return None


def _apply_required_native_dependency_config(config: dict) -> dict:
    for key, dependency in (
        ("cmake.define.EPCSAFT_ENABLE_CERES", "Ceres"),
        ("cmake.define.EPCSAFT_ENABLE_CPPAD", "CppAD"),
    ):
        value = _config_value(config, key)
        if value is not None and value.strip().upper() in {"0", "FALSE", "NO", "OFF"}:
            raise ValueError(f"{dependency} is required for native regression and derivative-capable package builds.")
        _set_config_default(config, key, "ON")
    return config


def _validate_ceres_dir(raw_path: str) -> Path:
    ceres_dir = Path(raw_path).expanduser().resolve()
    if not ceres_dir.is_dir():
        raise FileNotFoundError(f"EPCSAFT_PEP517_CERES_DIR does not exist or is not a directory: {ceres_dir}")
    if not any((ceres_dir / name).is_file() for name in ("CeresConfig.cmake", "ceres-config.cmake")):
        raise FileNotFoundError(
            "EPCSAFT_PEP517_CERES_DIR must point at the directory containing CeresConfig.cmake "
            f"or ceres-config.cmake: {ceres_dir}"
        )
    return ceres_dir


def _validate_ipopt_dir(raw_path: str) -> Path:
    ipopt_dir = Path(raw_path).expanduser().resolve()
    if not ipopt_dir.is_dir():
        raise FileNotFoundError(f"EPCSAFT_PEP517_IPOPT_DIR does not exist or is not a directory: {ipopt_dir}")
    if not any((ipopt_dir / name).is_file() for name in ("IpoptConfig.cmake", "ipopt-config.cmake")):
        raise FileNotFoundError(
            "EPCSAFT_PEP517_IPOPT_DIR must point at the directory containing IpoptConfig.cmake "
            f"or ipopt-config.cmake: {ipopt_dir}"
        )
    return ipopt_dir


def _validate_ipopt_root(raw_path: str) -> Path:
    ipopt_root = Path(raw_path).expanduser().resolve()
    if not ipopt_root.is_dir():
        raise FileNotFoundError(f"EPCSAFT_PEP517_IPOPT_ROOT does not exist or is not a directory: {ipopt_root}")
    include_dir = ipopt_root / "include"
    lib_dir = ipopt_root / "lib"
    if not include_dir.is_dir() or not lib_dir.is_dir():
        raise FileNotFoundError(
            "EPCSAFT_PEP517_IPOPT_ROOT must point at an Ipopt tree with include/ and lib/: "
            f"{ipopt_root}"
        )
    headers = (
        include_dir / "coin-or" / "IpIpoptApplication.hpp",
        include_dir / "coin" / "IpIpoptApplication.hpp",
        include_dir / "IpIpoptApplication.hpp",
    )
    libraries = (
        lib_dir / "ipopt.lib",
        lib_dir / "ipopt-3.lib",
        lib_dir / "libipopt.dll.a",
        lib_dir / "libipopt.a",
    )
    if not any(path.is_file() for path in headers):
        raise FileNotFoundError(f"EPCSAFT_PEP517_IPOPT_ROOT is missing Ipopt C++ headers: {ipopt_root}")
    if not any(path.is_file() for path in libraries):
        raise FileNotFoundError(f"EPCSAFT_PEP517_IPOPT_ROOT is missing an Ipopt link library: {ipopt_root}")
    return ipopt_root


def _apply_system_ceres_config(config: dict) -> dict:
    ceres_dir_env = os.environ.get("EPCSAFT_PEP517_CERES_DIR") or os.environ.get("Ceres_DIR")
    use_system_ceres = bool(ceres_dir_env) or _truthy_env("EPCSAFT_PEP517_USE_SYSTEM_CERES")
    if not use_system_ceres:
        return config

    _set_config_default(config, "cmake.define.EPCSAFT_ENABLE_CERES", "ON")
    _set_config_default(config, "cmake.define.EPCSAFT_USE_SYSTEM_CERES", "ON")
    if ceres_dir_env:
        _set_config_default(config, "cmake.define.Ceres_DIR", str(_validate_ceres_dir(ceres_dir_env)))
    return config


def _apply_system_ipopt_config(config: dict) -> dict:
    ipopt_dir_env = os.environ.get("EPCSAFT_PEP517_IPOPT_DIR") or os.environ.get("Ipopt_DIR")
    ipopt_root_env = os.environ.get("EPCSAFT_PEP517_IPOPT_ROOT") or os.environ.get("EPCSAFT_IPOPT_ROOT")
    if ipopt_dir_env and ipopt_root_env:
        raise ValueError("Use either EPCSAFT_PEP517_IPOPT_DIR or EPCSAFT_PEP517_IPOPT_ROOT, not both.")
    use_system_ipopt = (
        bool(ipopt_dir_env)
        or bool(ipopt_root_env)
        or _truthy_env("EPCSAFT_PEP517_ENABLE_IPOPT")
        or _truthy_env("EPCSAFT_PEP517_USE_SYSTEM_IPOPT")
    )
    if not use_system_ipopt:
        return config

    _set_config_default(config, "cmake.define.EPCSAFT_ENABLE_IPOPT", "ON")
    _set_config_default(config, "cmake.define.EPCSAFT_USE_SYSTEM_IPOPT", "ON")
    if ipopt_dir_env:
        _set_config_default(config, "cmake.define.Ipopt_DIR", str(_validate_ipopt_dir(ipopt_dir_env)))
    if ipopt_root_env:
        ipopt_root = _validate_ipopt_root(ipopt_root_env)
        _set_config_default(config, "cmake.define.EPCSAFT_IPOPT_ROOT", str(ipopt_root))
        bin_dir = ipopt_root / "bin"
        if bin_dir.is_dir():
            os.environ["PATH"] = str(bin_dir) + os.pathsep + os.environ.get("PATH", "")
            os.environ["EPCSAFT_RUNTIME_DLL_DIRS"] = (
                str(bin_dir) + os.pathsep + os.environ.get("EPCSAFT_RUNTIME_DLL_DIRS", "")
            )
    return config


def _apply_native_dependency_config(config: dict) -> dict:
    return _apply_system_ipopt_config(_apply_system_ceres_config(_apply_required_native_dependency_config(config)))


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
        return _apply_native_dependency_config(config)
    persistent = os.environ.get("EPCSAFT_PEP517_BUILD_DIR")
    if persistent:
        build_dir = Path(persistent).expanduser().resolve()
        build_dir.mkdir(parents=True, exist_ok=True)
    else:
        build_dir = Path(tempfile.mkdtemp(prefix="epcsaft-pep517-build-", dir=_external_temp_root())).resolve()
    config["build-dir"] = str(build_dir)
    return _apply_native_dependency_config(config)


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    return _scikit_build.build_wheel(
        wheel_directory,
        config_settings=_isolated_build_config(config_settings),
        metadata_directory=metadata_directory,
    )


def build_editable(wheel_directory, config_settings=None, metadata_directory=None):
    return _scikit_build.build_editable(
        wheel_directory,
        config_settings=_isolated_build_config(config_settings),
        metadata_directory=metadata_directory,
    )


build_sdist = _scikit_build.build_sdist
get_requires_for_build_editable = _scikit_build.get_requires_for_build_editable
get_requires_for_build_sdist = _scikit_build.get_requires_for_build_sdist
get_requires_for_build_wheel = _scikit_build.get_requires_for_build_wheel
prepare_metadata_for_build_editable = _scikit_build.prepare_metadata_for_build_editable
prepare_metadata_for_build_wheel = _scikit_build.prepare_metadata_for_build_wheel
