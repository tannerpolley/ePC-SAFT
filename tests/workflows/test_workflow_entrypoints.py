from __future__ import annotations

from pathlib import Path
import tomllib

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def test_bootstrap_scripts_use_normal_build_and_confidence_suite() -> None:
    for path in ("scripts/bootstrap_uv.ps1", "scripts/bootstrap_uv.sh"):
        content = _read(path)

        assert "uv python pin 3.12" in content
        assert "uv sync --no-install-project" in content
        assert "scripts/build_epcsaft.py --clean" not in content
        assert "scripts\\build_epcsaft.py --clean" not in content
        assert "run_pytest.py --confidence -q" in content
        assert "run_pytest.py tests/test_runtime.py -q" not in content
        assert "run_pytest.py tests\\test_runtime.py -q" not in content


def test_clean_scripts_announce_repair_only_scope() -> None:
    for path in ("scripts/clean_build.ps1", "scripts/clean_build.sh"):
        content = _read(path)

        assert "REPAIR-ONLY" in content
        assert "build/cache/native artifacts" in content


def test_codex_actions_name_standard_confidence_suite() -> None:
    data = tomllib.loads(_read(".codex/environments/environment.toml"))
    actions = {action["name"]: action["command"] for action in data["actions"]}

    assert "Test Standard Confidence Suite" in actions
    assert "run_pytest.py --confidence -q" in actions["Test Standard Confidence Suite"]
    assert "Test Confidence Suite (External Temp)" in actions
    assert actions["Fast Native Rebuild"] == "uv run python scripts/build_epcsaft.py --build-only --parallel 10"
    assert actions["Runtime Profile (Quick)"] == "uv run python run_pytest.py --profile -q"
    assert actions["Runtime Profile (Full)"] == "uv run python run_pytest.py --profile-full -q -s"
    assert actions["Build Distribution + Smoke Import"] == "uv run python scripts/build_dist.py"


def test_docs_make_confidence_suite_the_default_runtime_check() -> None:
    readme = _read("README.md")
    getting_started = _read("docs/pages/getting_started.rst")
    overview = _read("docs/pages/README.rst")
    docs_index = _read("docs/pages/index.rst")
    codex_workflows = _read("docs/pages/codex_workflows.rst")

    assert "default new-agent validation sequence" in readme
    assert "`--confidence` is the default runtime-confidence check" in readme
    assert "Codex workflow guide" in readme
    assert "default new-agent validation sequence" in getting_started
    assert "``--confidence`` is the default runtime-confidence check" in getting_started
    assert "uv run python run_pytest.py --confidence -q" in overview
    assert "run_pytest.py tests/test_runtime.py -q" not in overview
    assert "codex_workflows" in docs_index
    assert "native_debugging" in docs_index
    assert "native/equation debugging guide" in getting_started
    assert "Start every new Codex thread with this sequence" in codex_workflows
    assert "uv run python scripts/build_epcsaft.py --build-only --parallel 10" in codex_workflows
    assert "uv run python run_pytest.py --runtime -q" in codex_workflows
    assert "uv run python run_pytest.py --profile -q" in codex_workflows
    assert "uv run python run_pytest.py --profile-full -q -s" in codex_workflows
    assert "uv run python run_pytest.py --list-slices" in codex_workflows
    assert "EPCSAFT_PYTEST_TEMP_ROOT" in codex_workflows
    assert "reuse them inside hot loops" in codex_workflows
    assert "``--profile`` is the quick runtime-only profile" in codex_workflows
    assert "``--profile-full`` runs runtime, MIAC, and regression profiles" in codex_workflows
    assert "allow at least 120 seconds" in codex_workflows
    assert "uv run python scripts/build_dist.py" in codex_workflows
    assert "Do not use ``--clean`` for routine validation" in codex_workflows
