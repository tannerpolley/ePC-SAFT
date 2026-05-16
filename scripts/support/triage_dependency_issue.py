from __future__ import annotations

import argparse
import json
import re
import subprocess
from dataclasses import dataclass
from typing import Any

REPOSITORY = "tannerpolley/ePC-SAFT"

REQUIRED_REPORT_FIELDS = (
    "Downstream repo or path",
    "Task goal",
    "Failing command",
    "Error or bad result",
    "Minimal reproducer",
    "Imported epcsaft path",
    "Expected behavior",
    "Actual behavior",
    "Downstream validation command after fix",
    "Temporary workaround",
)
GATING_FIELDS = (
    "Minimal reproducer",
    "Imported epcsaft path",
    "Downstream validation command after fix",
)
AREA_LABELS = (
    "python-api",
    "native",
    "solver",
    "packaging",
    "docs",
    "validation",
    "regression",
    "equilibrium",
)
BLOCKING_LABELS = {"in-progress", "blocked-downstream", "downstream-validated"}

KEYWORDS = {
    "packaging": ("install", "wheel", "build", "pip", "uv", "_core", "import path", "editable"),
    "native": ("native", "c++", "pybind", "_core", "segfault", "density", "helmholtz"),
    "solver": ("solver", "converge", "residual", "flash", "bubble", "stability"),
    "equilibrium": ("equilibrium", "lle", "vle", "phase", "speciation", "bubble"),
    "regression": ("fit", "regression", "objective", "parameter", "least" + "_squares"),
    "validation": ("validation", "test", "benchmark", "expected", "mismatch"),
    "docs": ("docs", "documentation", "readme", "ambiguous"),
    "python-api": ("api", "public", "method", "argument", "keyword", "attribute"),
}


@dataclass(frozen=True)
class IssueTriage:
    issue: dict[str, Any]
    fields: dict[str, str]
    missing_fields: list[str]
    labels: list[str]
    classification: str
    recommended_commands: list[str]
    ready: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "issue": self.issue,
            "fields": self.fields,
            "missing_fields": self.missing_fields,
            "labels": self.labels,
            "classification": self.classification,
            "recommended_commands": self.recommended_commands,
            "ready": self.ready,
        }


def issue_number_from_token(token: str) -> str:
    value = str(token).strip()
    match = re.search(r"/issues/(\d+)(?:\D.*)?$", value)
    if match:
        return match.group(1)
    if re.fullmatch(r"\d+", value):
        return value
    raise SystemExit(f"Could not parse GitHub issue number from: {token}")


def fetch_issue(issue: str) -> dict[str, Any]:
    issue_number = issue_number_from_token(issue)
    cmd = [
        "gh",
        "issue",
        "view",
        issue_number,
        "--repo",
        REPOSITORY,
        "--json",
        "number,title,body,labels,url",
    ]
    completed = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "gh issue view failed"
        raise SystemExit(message)
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"gh issue view returned invalid JSON: {exc}") from exc


def parse_issue_fields(body: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    current: str | None = None
    lines: list[str] = []
    for line in str(body or "").splitlines():
        heading = re.match(r"^###\s+(.+?)\s*$", line)
        if heading:
            if current is not None:
                fields[current] = _clean_field_value("\n".join(lines))
            current = heading.group(1).strip()
            lines = []
            continue
        if current is not None:
            lines.append(line)
    if current is not None:
        fields[current] = _clean_field_value("\n".join(lines))
    return fields


def triage_issue(issue_payload: dict[str, Any]) -> IssueTriage:
    body = str(issue_payload.get("body") or "")
    fields = parse_issue_fields(body)
    labels = _labels_from_issue(issue_payload)
    missing = [field for field in REQUIRED_REPORT_FIELDS if not _has_value(fields.get(field))]
    ready = not any(field in missing for field in GATING_FIELDS)
    classification = classify_issue(labels, body)
    recommended = recommended_commands(issue_payload, missing, labels, classification)
    issue = {
        "number": issue_payload.get("number"),
        "title": issue_payload.get("title"),
        "url": issue_payload.get("url"),
    }
    return IssueTriage(
        issue=issue,
        fields=fields,
        missing_fields=missing,
        labels=labels,
        classification=classification,
        recommended_commands=recommended,
        ready=ready,
    )


def classify_issue(labels: list[str], body: str) -> str:
    for label in labels:
        if label in AREA_LABELS:
            return label
    haystack = body.lower()
    scores = {area: sum(1 for keyword in keywords if keyword in haystack) for area, keywords in KEYWORDS.items()}
    best = max(scores, key=lambda area: scores[area])
    return best if scores[best] > 0 else "python-api"


def recommended_commands(
    issue_payload: dict[str, Any],
    missing_fields: list[str],
    labels: list[str],
    classification: str,
) -> list[str]:
    number = issue_payload.get("number", "<issue>")
    commands: list[str] = []
    if BLOCKING_LABELS.intersection(labels):
        commands.append(f"gh issue view {number} --repo {REPOSITORY}")
        return commands
    if any(field in missing_fields for field in GATING_FIELDS):
        commands.append(f"gh issue edit {number} --repo {REPOSITORY} --add-label needs-repro")
        return commands
    commands.append(f"uv run python scripts/support/triage_dependency_issue.py --issue {number} --json")
    commands.append("uv run python run_pytest.py <focused-test-targets> -q")
    if classification in {"native", "packaging"}:
        commands.append("uv run python scripts/dev/build_epcsaft.py")
        commands.append("uv run python scripts/dev/doctor.py")
    commands.append("uv run python scripts/dev/validate_project.py quick")
    return commands


def render_text(triage: IssueTriage) -> str:
    lines = [
        f"Issue: #{triage.issue.get('number')} {triage.issue.get('title')}",
        f"URL: {triage.issue.get('url')}",
        f"Labels: {', '.join(triage.labels) if triage.labels else '<none>'}",
        f"Classification: {triage.classification}",
        f"Ready: {triage.ready}",
    ]
    if triage.missing_fields:
        lines.append("Missing fields:")
        lines.extend(f"- {field}" for field in triage.missing_fields)
    lines.append("Recommended commands:")
    lines.extend(f"- {command}" for command in triage.recommended_commands)
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Read-only triage for downstream epcsaft dependency issues.")
    parser.add_argument("--issue", required=True, help="GitHub issue number or URL.")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable triage JSON.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    triage = triage_issue(fetch_issue(args.issue))
    if args.json:
        print(json.dumps(triage.to_dict(), indent=2, sort_keys=True))
    else:
        print(render_text(triage))
    return 0 if triage.ready else 1


def _labels_from_issue(issue: dict[str, Any]) -> list[str]:
    labels = []
    for item in issue.get("labels") or []:
        if isinstance(item, dict):
            name = item.get("name")
        else:
            name = item
        if name:
            labels.append(str(name))
    return labels


def _clean_field_value(value: str) -> str:
    return re.sub(r"\n{3,}", "\n\n", value).strip()


def _has_value(value: str | None) -> bool:
    if value is None:
        return False
    stripped = value.strip()
    return bool(stripped) and stripped.lower() not in {"none", "n/a", "na", "_no response_"}


if __name__ == "__main__":
    raise SystemExit(main())
