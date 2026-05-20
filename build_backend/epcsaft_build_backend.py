"""PEP 517 backend wrapper for sandbox-safe Windows package builds."""

from __future__ import annotations

import errno
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from scikit_build_core import build as _scikit_build

try:
    from native_dependency_policy import (
        ipopt_root_prefers_msvc,
        resolve_default_windows_ipopt_sdk_root,
        validate_ceres_dir,
        validate_ipopt_dir,
        validate_ipopt_root,
    )
except ModuleNotFoundError:  # pragma: no cover - import-by-file tests
    _POLICY_PATH = Path(__file__).with_name("native_dependency_policy.py")
    _POLICY_SPEC = importlib.util.spec_from_file_location("epcsaft_native_dependency_policy", _POLICY_PATH)
    if _POLICY_SPEC is None or _POLICY_SPEC.loader is None:
        raise
    _POLICY = importlib.util.module_from_spec(_POLICY_SPEC)
    _POLICY_SPEC.loader.exec_module(_POLICY)
    ipopt_root_prefers_msvc = _POLICY.ipopt_root_prefers_msvc
    resolve_default_windows_ipopt_sdk_root = _POLICY.resolve_default_windows_ipopt_sdk_root
    validate_ceres_dir = _POLICY.validate_ceres_dir
    validate_ipopt_dir = _POLICY.validate_ipopt_dir
    validate_ipopt_root = _POLICY.validate_ipopt_root


def _find_vsdevcmd() -> Path:
    explicit = os.environ.get("EPCSAFT_VSDEVCMD")
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if path.is_file():
            return path
    vswhere = Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")) / (
        "Microsoft Visual Studio/Installer/vswhere.exe"
    )
    if vswhere.is_file():
        output = subprocess.check_output(
            [
                str(vswhere),
                "-latest",
                "-products",
                "*",
                "-requires",
                "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                "-property",
                "installationPath",
            ],
            text=True,
        ).strip()
        if output:
            candidate = Path(output) / "Common7" / "Tools" / "VsDevCmd.bat"
            if candidate.is_file():
                return candidate
    default_path = Path(r"C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat")
    if default_path.is_file():
        return default_path
    raise FileNotFoundError("MSVC Ipopt SDK requested, but VsDevCmd.bat was not found.")


def _load_msvc_env_for_ipopt_root(ipopt_root: Path) -> None:
    if not ipopt_root_prefers_msvc(ipopt_root) or shutil.which("cl"):
        return
    if os.environ.get("CMAKE_GENERATOR") == "MinGW Makefiles":
        raise RuntimeError("The MSVC Ipopt SDK cannot be used with the MinGW Makefiles generator.")
    vsdevcmd = _find_vsdevcmd()
    command = f'call "{vsdevcmd}" -arch=x64 -host_arch=x64 >nul && set'
    output = subprocess.check_output(command, text=True, shell=True)
    for line in output.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        os.environ[key] = value


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


def _apply_system_ceres_config(config: dict) -> dict:
    ceres_dir_env = os.environ.get("EPCSAFT_PEP517_CERES_DIR") or os.environ.get("Ceres_DIR")
    use_system_ceres = bool(ceres_dir_env) or _truthy_env("EPCSAFT_PEP517_USE_SYSTEM_CERES")
    if not use_system_ceres:
        return config

    _set_config_default(config, "cmake.define.EPCSAFT_ENABLE_CERES", "ON")
    _set_config_default(config, "cmake.define.EPCSAFT_USE_SYSTEM_CERES", "ON")
    if ceres_dir_env:
        _set_config_default(config, "cmake.define.Ceres_DIR", str(validate_ceres_dir(ceres_dir_env)))
    return config


def _apply_system_ipopt_config(config: dict) -> dict:
    ipopt_dir_env = os.environ.get("EPCSAFT_PEP517_IPOPT_DIR") or os.environ.get("Ipopt_DIR")
    ipopt_root_env = os.environ.get("EPCSAFT_PEP517_IPOPT_ROOT") or os.environ.get("EPCSAFT_IPOPT_ROOT")
    if not ipopt_dir_env and not ipopt_root_env:
        default_ipopt_root = resolve_default_windows_ipopt_sdk_root()
        if default_ipopt_root is not None:
            ipopt_root_env = str(default_ipopt_root)
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
        _set_config_default(config, "cmake.define.Ipopt_DIR", str(validate_ipopt_dir(ipopt_dir_env)))
    if ipopt_root_env:
        ipopt_root = validate_ipopt_root(ipopt_root_env)
        _load_msvc_env_for_ipopt_root(ipopt_root)
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
