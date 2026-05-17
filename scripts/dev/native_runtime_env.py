from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import MutableMapping

REPO_ROOT = Path(__file__).resolve().parents[2]
DEV_BUILD_CACHE = REPO_ROOT / "build" / "dev" / "CMakeCache.txt"


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


def resolve_ipopt_root(
    *,
    cache_path: Path = DEV_BUILD_CACHE,
    explicit_root: Path | str | None = None,
    env: MutableMapping[str, str] | None = None,
) -> Path | None:
    runtime_env = os.environ if env is None else env
    raw = (
        str(explicit_root)
        if explicit_root is not None
        else cmake_cache_value("EPCSAFT_IPOPT_ROOT", cache_path)
        or runtime_env.get("EPCSAFT_IPOPT_ROOT")
        or runtime_env.get("EPCSAFT_PEP517_IPOPT_ROOT")
    )
    if raw is None:
        return None
    raw = raw.strip()
    if not raw or raw == "<unconfigured>":
        return None
    return Path(raw).expanduser().resolve()


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
        _prepend_unique_path(runtime_env, "PATH", runtime_dir)
        _prepend_unique_path(runtime_env, "EPCSAFT_RUNTIME_DLL_DIRS", runtime_dir)
        applied = True
    return NativeRuntimeEnv(
        ipopt_configured=configured,
        ipopt_root=root,
        ipopt_runtime_dir=runtime_dir,
        applied=applied,
    )


def apply_to_current_process(*, cache_path: Path = DEV_BUILD_CACHE) -> NativeRuntimeEnv:
    return apply_native_runtime_env(os.environ, cache_path=cache_path)
