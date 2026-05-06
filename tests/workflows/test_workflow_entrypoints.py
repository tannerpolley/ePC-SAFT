from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def test_bootstrap_scripts_use_normal_build_and_fast_suite() -> None:
    for path in ("scripts/bootstrap_uv.ps1", "scripts/bootstrap_uv.sh"):
        content = _read(path)

        assert "uv python pin 3.12" in content
        assert "uv sync --no-install-project" in content
        assert "scripts/build_epcsaft.py --clean" not in content
        assert "scripts\\build_epcsaft.py --clean" not in content
        assert "scripts\\validate_project.py quick" in content or "scripts/validate_project.py quick" in content
        assert "run_pytest.py tests/test_runtime.py -q" not in content
        assert "run_pytest.py tests\\test_runtime.py -q" not in content


def test_clean_scripts_announce_repair_only_scope() -> None:
    for path in ("scripts/clean_build.ps1", "scripts/clean_build.sh"):
        content = _read(path)

        assert "REPAIR-ONLY" in content
        assert "build/cache/native artifacts" in content


def test_docs_make_confidence_suite_the_default_runtime_check() -> None:
    readme = _read("README.md")
    getting_started = _read("docs/pages/getting_started.rst")
    overview = _read("docs/pages/README.rst")
    docs_index = _read("docs/pages/index.rst")
    development_workflows = _read("docs/pages/development_workflows.rst")

    assert "default source-checkout validation sequence" in readme
    assert "`run_pytest.py -q` is the default fast contract suite" in readme
    assert "development workflow guide" in readme
    assert "default source-checkout validation sequence" in getting_started
    assert "``run_pytest.py -q`` is the default fast contract suite" in getting_started
    assert "uv run python run_pytest.py --confidence -q" in overview
    assert "run_pytest.py tests/test_runtime.py -q" not in overview
    assert "development_workflows" in docs_index
    assert "native_debugging" in docs_index
    assert "native/equation debugging guide" in getting_started
    assert "Start every fresh source checkout with this sequence" in development_workflows
    assert "uv run python scripts/build_epcsaft.py --build-only --parallel 10" in development_workflows
    assert "uv run python run_pytest.py --runtime -q" in development_workflows
    assert "uv run python run_pytest.py --profile -q" in development_workflows
    assert "uv run python run_pytest.py --profile-full -q -s" in development_workflows
    assert "uv run python run_pytest.py --list-slices" in development_workflows
    assert "EPCSAFT_PYTEST_TEMP_ROOT" in development_workflows
    assert "reuse them inside hot loops" in development_workflows
    assert "``--profile`` is the quick runtime-only profile" in development_workflows
    assert "``--profile-full`` runs runtime, MIAC, and regression profiles" in development_workflows
    assert "allow at least 120 seconds" in development_workflows
    assert "uv run python scripts/build_dist.py" in development_workflows
    assert "Do not use ``--clean`` for routine validation" in development_workflows
