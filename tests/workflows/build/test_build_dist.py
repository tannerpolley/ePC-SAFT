from __future__ import annotations

import zipfile

import pytest

from scripts.dev import build_dist


def _wheel(path, names: list[str]) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        for name in names:
            archive.writestr(name, "")


def test_dist_wheel_audit_accepts_package_only_wheel(tmp_path) -> None:
    wheel = tmp_path / "epcsaft-0-py3-none-any.whl"
    _wheel(wheel, ["epcsaft/__init__.py", "epcsaft/_core.cp313-win_amd64.pyd"])

    build_dist._assert_wheel_has_no_solver_dev_artifacts(wheel)


def test_dist_wheel_audit_rejects_ceres_development_artifacts(tmp_path) -> None:
    wheel = tmp_path / "epcsaft-0-py3-none-any.whl"
    _wheel(
        wheel,
        [
            "include/ceres/solver.h",
            "lib/cmake/Ceres/CeresConfig.cmake",
            "include/ceres/" + "numeric" + "_diff_cost_function.h",
        ],
    )

    with pytest.raises(RuntimeError, match="Ceres development artifacts"):
        build_dist._assert_wheel_has_no_solver_dev_artifacts(wheel)


def test_dist_build_env_sets_conservative_parallel_level(monkeypatch) -> None:
    monkeypatch.delenv("PYTHONPATH", raising=False)

    env = build_dist._env("1")

    assert env["CMAKE_BUILD_PARALLEL_LEVEL"] == "1"
    assert env["EPCSAFT_SANDBOX_SAFE_TEMPFILE"] == "1"
