from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def test_bootstrap_scripts_use_normal_build_and_fast_suite() -> None:
    for path in ("scripts/bootstrap_uv.ps1", "scripts/bootstrap_uv.sh"):
        content = _read(path)

        assert "uv python pin 3.13" in content
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
    release_installation = _read("docs/pages/release_installation.rst")
    docs_index = _read("docs/pages/index.rst")
    development_workflows = _read("docs/pages/development_workflows.rst")

    assert "README intentionally stays focused on package users" in readme
    assert "uv run python scripts\\validate_project.py quick" not in readme
    assert "uv run python run_pytest.py --confidence -q" not in readme
    assert "The current public release is also available from GitHub" in readme
    assert "python -m pip install epcsaft" in readme
    assert "python -m pip install -e ." in readme
    assert "README intentionally stays focused on package users" in readme
    assert "Editable source install" in release_installation
    assert "python -m pip install -e ." in release_installation
    assert "Source and editable installs build a native C++ extension" in getting_started
    assert "default source-checkout validation sequence" not in getting_started
    assert "``run_pytest.py -q`` is the default fast contract suite" not in getting_started
    assert "Current release: ``1.5.2``" in overview
    assert "python -m pip install -e ." in overview
    assert "uv run python run_pytest.py --confidence -q" not in overview
    assert "run_pytest.py tests/test_runtime.py -q" not in overview
    assert "release_installation" in docs_index
    assert "development_workflows" in docs_index
    assert "native_debugging" in docs_index
    assert "publishing" in docs_index
    assert "native/equation debugging guide" not in getting_started
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


def test_github_default_smoke_uses_downstream_path_install_not_wheel_build() -> None:
    workflow = _read(".github/workflows/wheels.yml")
    old_python = "uv python install " + "3." + "12"
    old_wheel_command = "uv build --" + "wheel"
    old_wheel_step = "Build Windows " + "wheel"

    assert "uv python install 3.13" in workflow
    assert old_python not in workflow
    assert '$repoUrl = "file:///"' in workflow
    assert "epcsaft @ $repoUrl" in workflow
    assert "UV_CACHE_DIR" in workflow
    assert "EPCSAFT_PEP517_BUILD_DIR" in workflow
    assert "uv sync --python 3.13" in workflow
    assert "uv run --no-sync python" in workflow
    assert old_wheel_command not in workflow
    assert old_wheel_step not in workflow


def test_github_full_packaging_remains_manual_only() -> None:
    workflow = _read(".github/workflows/wheels.yml")
    old_cibw = 'CIBW_BUILD: "cp' + '312-*"'

    assert "workflow_dispatch:" in workflow
    assert "if: ${{ github.event_name == 'workflow_dispatch' && inputs.full_wheel_matrix }}" in workflow
    assert 'CIBW_BUILD: "cp313-*"' in workflow
    assert old_cibw not in workflow


def test_github_default_events_do_not_run_duplicate_heavy_smokes() -> None:
    workflow = _read(".github/workflows/wheels.yml")

    assert "fast-pr-smoke:" in workflow
    assert "windows-install-smoke:" in workflow
    assert "if: ${{ github.event_name == 'pull_request' }}" in workflow
    assert "if: ${{ github.event_name != 'pull_request' }}" in workflow
    assert workflow.count("name: windows install smoke") == 1
    assert workflow.count("name: fast workflow smoke") == 1


def test_pypi_publish_workflow_uses_trusted_publishing() -> None:
    workflow = _read(".github/workflows/publish-pypi.yml")
    publishing_docs = _read("docs/pages/publishing.rst")

    for token in (
        "name: publish-to-pypi",
        "types: [published]",
        "workflow_dispatch:",
        "id-token: write",
        "environment:",
        "name: pypi",
        "pypa/gh-action-pypi-publish@release/v1",
        "actions/download-artifact@v8.0.1",
        "merge-multiple: true",
        'CIBW_BUILD: "cp313-*"',
        "uv build --sdist",
    ):
        assert token in workflow
    assert "password:" not in workflow
    assert "username:" not in workflow
    assert "PYPI_API_TOKEN" not in workflow
    assert "pp*" not in workflow
    assert "Workflow filename: ``publish-pypi.yml``" in publishing_docs
    assert "Environment name: ``pypi``" in publishing_docs
