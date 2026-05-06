from __future__ import annotations

import json
import subprocess
from pathlib import Path

import yaml

from scripts import triage_dependency_issue as triage

REPO_ROOT = Path(__file__).resolve().parents[2]
ISSUE_FORMS = (
    REPO_ROOT / ".github" / "ISSUE_TEMPLATE" / "downstream_dependency_bug.yml",
    REPO_ROOT / ".github" / "ISSUE_TEMPLATE" / "upstream_package_request.yml",
)
REQUIRED_LABELS = {
    "downstream-bug",
    "needs-repro",
    "agent-ready",
    "in-progress",
    "upstream-fix-ready",
    "downstream-validated",
    "blocked-downstream",
    "python-api",
    "native",
    "solver",
    "packaging",
    "docs",
    "validation",
    "regression",
    "equilibrium",
}


def _body(**overrides: str) -> str:
    values = {
        "Downstream repo or path": "C:\\Users\\Tanner\\Documents\\git\\MEA-Thermodynamics",
        "Task goal": "Run electrolyte VLE against epcsaft.",
        "Failing command": "uv run python analysis\\run_case.py",
        "Error or bad result": "SolutionError: residual did not converge",
        "Minimal reproducer": "import epcsaft\nprint(epcsaft.__file__)",
        "Imported epcsaft path": "C:\\Users\\Tanner\\Documents\\git\\ePC-SAFT\\src\\epcsaft\\__init__.py",
        "Expected behavior": "The public API should solve the case.",
        "Actual behavior": "The solver fails before returning diagnostics.",
        "Downstream validation command after fix": "uv run python analysis\\run_case.py",
        "Temporary workaround": "None.",
    }
    values.update(overrides)
    return "\n\n".join(f"### {key}\n\n{value}" for key, value in values.items())


def _issue_payload(body: str, labels: list[str] | None = None) -> dict[str, object]:
    return {
        "number": 12,
        "title": "Downstream failure: MEA hits solver diagnostics",
        "url": "https://github.com/tannerpolley/ePC-SAFT/issues/12",
        "body": body,
        "labels": [{"name": label} for label in (labels or ["downstream-bug", "agent-ready"])],
    }


def test_issue_forms_are_valid_yaml_and_require_downstream_report_fields() -> None:
    required_ids = {
        "downstream-repo-or-path",
        "task-goal",
        "failing-command",
        "error-or-bad-result",
        "minimal-reproducer",
        "imported-epcsaft-path",
        "expected-behavior",
        "actual-behavior",
        "downstream-validation-command-after-fix",
        "temporary-workaround",
    }

    for path in ISSUE_FORMS:
        form = yaml.safe_load(path.read_text(encoding="utf-8"))
        body = form["body"]
        fields = {item.get("id"): item for item in body if "id" in item}

        assert required_ids.issubset(fields), path
        for field_id in required_ids:
            assert fields[field_id]["validations"]["required"] is True


def test_issue_form_labels_match_planned_taxonomy() -> None:
    label_text = "\n".join(path.read_text(encoding="utf-8") for path in ISSUE_FORMS)

    downstream_labels = yaml.safe_load(ISSUE_FORMS[0].read_text(encoding="utf-8"))["labels"]
    request_labels = yaml.safe_load(ISSUE_FORMS[1].read_text(encoding="utf-8"))["labels"]

    assert {"downstream-bug", "agent-ready"}.issubset(downstream_labels)
    assert {"enhancement", "agent-ready"}.issubset(request_labels)
    assert REQUIRED_LABELS.issuperset(triage.AREA_LABELS)
    assert {"downstream-bug", "needs-repro"}.issubset(REQUIRED_LABELS)
    assert "agent-ready" in REQUIRED_LABELS
    assert "blocked-downstream" in REQUIRED_LABELS
    assert "docs" in label_text


def test_complete_downstream_issue_returns_ready_triage() -> None:
    result = triage.triage_issue(_issue_payload(_body(), labels=["downstream-bug", "agent-ready", "solver"]))

    assert result.ready is True
    assert result.missing_fields == []
    assert result.classification == "solver"
    assert "uv run python scripts/validate_project.py quick" in result.recommended_commands
    assert result.fields["Imported epcsaft path"].endswith("epcsaft\\__init__.py")


def test_missing_minimal_reproducer_returns_needs_repro_recommendation() -> None:
    result = triage.triage_issue(_issue_payload(_body(**{"Minimal reproducer": ""}), labels=["downstream-bug"]))

    assert result.ready is False
    assert "Minimal reproducer" in result.missing_fields
    assert result.recommended_commands == ["gh issue edit 12 --repo tannerpolley/ePC-SAFT --add-label needs-repro"]


def test_issue_number_and_url_resolve_to_same_gh_issue_view_call(monkeypatch) -> None:
    calls: list[list[str]] = []

    def fake_run(cmd, check, capture_output, text):
        calls.append(cmd)
        payload = json.dumps(_issue_payload(_body()))
        return subprocess.CompletedProcess(cmd, 0, stdout=payload, stderr="")

    monkeypatch.setattr(triage.subprocess, "run", fake_run)

    assert triage.fetch_issue("12")["number"] == 12
    assert triage.fetch_issue("https://github.com/tannerpolley/ePC-SAFT/issues/12")["number"] == 12
    assert calls[0] == calls[1]
    assert calls[0][:4] == ["gh", "issue", "view", "12"]


def test_json_output_is_machine_readable(monkeypatch, capsys) -> None:
    def fake_fetch(issue: str) -> dict[str, object]:
        assert issue == "12"
        return _issue_payload(_body(), labels=["downstream-bug", "agent-ready", "packaging"])

    monkeypatch.setattr(triage, "fetch_issue", fake_fetch)

    exit_code = triage.main(["--issue", "12", "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["issue"]["number"] == 12
    assert payload["missing_fields"] == []
    assert payload["labels"] == ["downstream-bug", "agent-ready", "packaging"]
    assert payload["classification"] == "packaging"
    assert "recommended_commands" in payload
