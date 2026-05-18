from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import MutableMapping

REPO_ROOT = Path(__file__).resolve().parents[2]
DEV_BUILD_CACHE = REPO_ROOT / "build" / "dev" / "CMakeCache.txt"
DEFAULT_WINDOWS_IPOPT_SDK_RELATIVE = Path("Documents") / "deps" / "ipopt-msvc"


@dataclass(frozen=True)
class NativeRuntimeEnv:
    ipopt_configured: bool
    ipopt_root: Path | None
    ipopt_runtime_dir: Path | None
    applied: bool


def cmake_cache_value(name: str, cache_path: Path = DEV_BUILD_CACHE) -> str | None:
    if not cache_path.exists():
        return None
    prefix = f"{name}:"
    for line in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix):
            return line.split("=", 1)[1].strip()
    return None


def cmake_enabled(value: str | None) -> bool:
    return str(value or "").strip().upper() in {"1", "ON", "TRUE", "YES"}


def default_windows_ipopt_sdk_root(home: Path | None = None) -> Path:
    return (Path.home() if home is None else home).expanduser() / DEFAULT_WINDOWS_IPOPT_SDK_RELATIVE


def resolve_default_windows_ipopt_sdk_root(
    *,
    home: Path | None = None,
    platform_name: str | None = None,
) -> Path | None:
    if (os.name if platform_name is None else platform_name) != "nt":
        return None
    root = default_windows_ipopt_sdk_root(home)
    return root.resolve() if root.is_dir() else None


def resolve_ipopt_root(
    *,
    cache_path: Path = DEV_BUILD_CACHE,
    explicit_root: Path | str | None = None,
    env: MutableMapping[str, str] | None = None,
) -> Path | None:
    runtime_env = os.environ if env is None else env
    raw = str(explicit_root) if explicit_root is not None else None
    if raw is None:
        raw = (
            cmake_cache_value("EPCSAFT_IPOPT_ROOT", cache_path)
            or runtime_env.get("EPCSAFT_IPOPT_ROOT")
            or runtime_env.get("EPCSAFT_PEP517_IPOPT_ROOT")
        )
    if raw is None:
        if cmake_cache_value("Ipopt_DIR", cache_path):
            return None
        return resolve_default_windows_ipopt_sdk_root()
    raw = raw.strip()
    if not raw or raw == "<unconfigured>":
        if cmake_cache_value("Ipopt_DIR", cache_path):
            return None
        return resolve_default_windows_ipopt_sdk_root()
    return Path(raw).expanduser().resolve()


def resolve_ipopt_root_for_build(
    raw_path: Path | str | None,
    *,
    enable_ipopt: bool,
    ipopt_dir: Path | str | None = None,
    default_root: Path | str | None = None,
    label: str = "Ipopt root",
) -> Path | None:
    if raw_path is None and enable_ipopt and ipopt_dir is None:
        raw_path = default_root if default_root is not None else resolve_default_windows_ipopt_sdk_root()
    if raw_path is None:
        return None
    path = Path(raw_path).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f"{label} does not exist or is not a directory: {path}")
    return path


def prepend_ipopt_runtime_env(env: MutableMapping[str, str], ipopt_root: Path | str | None) -> Path | None:
    root = Path(ipopt_root).expanduser().resolve() if ipopt_root is not None else None
    runtime_dir = ipopt_runtime_bin(root)
    if runtime_dir is None:
        return None
    _prepend_unique_path(env, "PATH", runtime_dir)
    _prepend_unique_path(env, "EPCSAFT_RUNTIME_DLL_DIRS", runtime_dir)
    return runtime_dir


def ipopt_runtime_bin(ipopt_root: Path | None) -> Path | None:
    if ipopt_root is None:
        return None
    bin_dir = ipopt_root / "bin"
    return bin_dir if bin_dir.is_dir() else None


def _prepend_unique_path(env: MutableMapping[str, str], name: str, path: Path) -> None:
    entry = str(path.resolve())
    current = env.get(name, "")
    existing = [part for part in current.split(os.pathsep) if part]
    normalized_entry = os.path.normcase(os.path.normpath(entry))
    kept = [part for part in existing if os.path.normcase(os.path.normpath(part)) != normalized_entry]
    env[name] = os.pathsep.join([entry, *kept])


def apply_native_runtime_env(
    env: MutableMapping[str, str] | None = None,
    *,
    cache_path: Path = DEV_BUILD_CACHE,
    ipopt_root: Path | str | None = None,
    ipopt_enabled: bool | None = None,
) -> NativeRuntimeEnv:
    runtime_env = os.environ if env is None else env
    configured = cmake_enabled(cmake_cache_value("EPCSAFT_ENABLE_IPOPT", cache_path))
    if ipopt_enabled is not None:
        configured = bool(ipopt_enabled)

    root = resolve_ipopt_root(cache_path=cache_path, explicit_root=ipopt_root, env=runtime_env)
    runtime_dir = ipopt_runtime_bin(root)
    applied = False
    if configured and runtime_dir is not None:
        applied = prepend_ipopt_runtime_env(runtime_env, root) is not None
    return NativeRuntimeEnv(
        ipopt_configured=configured,
        ipopt_root=root,
        ipopt_runtime_dir=runtime_dir,
        applied=applied,
    )


def apply_to_current_process(*, cache_path: Path = DEV_BUILD_CACHE) -> NativeRuntimeEnv:
    return apply_native_runtime_env(os.environ, cache_path=cache_path)
