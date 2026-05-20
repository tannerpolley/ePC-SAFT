from __future__ import annotations

import argparse
from collections import defaultdict
import difflib
import json
import re
import shutil
import sys
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

TEX_PATH = REPO_ROOT / "docs" / "latex" / "algorithms.tex"
MARKDOWN_PATH = REPO_ROOT / "docs" / "algorithms.md"
REGISTRY_PATH = REPO_ROOT / "docs" / "algorithms_registry.yaml"
LEGACY_RENDER_DIR = REPO_ROOT / "docs" / "rendered_math" / "algorithms"
LEGACY_RENDER_ROOT = REPO_ROOT / "docs" / "rendered_math"
SCAN_ROOTS = (
    REPO_ROOT / "src" / "epcsaft",
    REPO_ROOT / "tests",
)

SECTION_RE = re.compile(r"\\section\{(.+?)\}")
SUBSECTION_RE = re.compile(r"\\subsection\{(.+?)\}")
SUBSUBSECTION_RE = re.compile(r"\\subsubsection\{(.+?)\}")
ALGID_RE = re.compile(r"%\s*AlgID:\s*([A-Za-z0-9_]+)")
META_RE = re.compile(r"%\s*([A-Za-z][A-Za-z /-]+):\s*(.*)")
CODE_ALGID_RE = re.compile(r"(?://|#)\s*AlgID:\s*([A-Za-z0-9_]+)")
DISPLAY_MATH_RE = re.compile(r"\\\[(.+?)\\\]", re.DOTALL)
INLINE_MATH_RE = re.compile(r"\\\((.+?)\\\)", re.DOTALL)
TEXTTT_RE = re.compile(r"\\texttt\{([^{}]+)\}")
SOURCE_SUFFIXES = {".cpp", ".h", ".hpp", ".py"}
DOCUMENTATION_ONLY_STATUSES = {
    "documentation-only",
    "documentation_only",
    "docs-only",
    "docs_only",
    "reference-only",
    "reference_only",
    "planned",
    "planned-only",
    "planned_only",
    "roadmap-only",
    "roadmap_only",
}


def metadata_key(label: str) -> str:
    return label.strip().lower().replace("-", "_").replace("/", "_").replace(" ", "_")


def _section_match(line: str) -> re.Match[str] | None:
    return SECTION_RE.search(line) or SUBSECTION_RE.search(line) or SUBSUBSECTION_RE.search(line)


def parse_algorithms(tex_path: Path) -> list[dict]:
    lines = tex_path.read_text(encoding="utf-8").splitlines()
    entries: list[dict] = []
    section = ""
    subsection = ""
    subsubsection = ""
    seen_algids: set[str] = set()
    i = 0
    while i < len(lines):
        line = lines[i]
        section_match = SECTION_RE.search(line)
        if section_match:
            section = section_match.group(1).strip()
            subsection = ""
            subsubsection = ""
            i += 1
            continue

        subsection_match = SUBSECTION_RE.search(line)
        if subsection_match:
            subsection = subsection_match.group(1).strip()
            subsubsection = ""
            i += 1
            continue

        subsubsection_match = SUBSUBSECTION_RE.search(line)
        if subsubsection_match:
            subsubsection = subsubsection_match.group(1).strip()
            i += 1
            continue

        algid_match = ALGID_RE.match(line.strip())
        if not algid_match:
            i += 1
            continue

        algid = algid_match.group(1)
        if algid in seen_algids:
            raise ValueError(f"Duplicate AlgID in {tex_path}: {algid}")
        seen_algids.add(algid)

        entry: dict[str, object] = {
            "algid": algid,
            "section": section,
            "subsection": subsection,
            "subsubsection": subsubsection,
            "tex_file": tex_path.relative_to(REPO_ROOT).as_posix(),
            "tex_line": i + 1,
        }

        i += 1
        while i < len(lines):
            stripped = lines[i].strip()
            meta_match = META_RE.match(stripped)
            if not meta_match:
                break
            entry[metadata_key(meta_match.group(1))] = meta_match.group(2).strip()
            i += 1

        body_lines: list[str] = []
        while i < len(lines):
            stripped = lines[i].strip()
            if stripped == r"\end{document}" or ALGID_RE.match(stripped) or _section_match(lines[i]):
                break
            body_lines.append(lines[i].rstrip())
            i += 1

        entry["body"] = "\n".join(body_lines).strip()
        entry["owner_refs"] = []
        entries.append(entry)

    return entries


def require_algorithms_tex(tex_path: Path = TEX_PATH) -> None:
    """Exit with an actionable message when the tracked LaTeX source is absent."""
    if tex_path.exists():
        return
    try:
        rel_path = tex_path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        rel_path = tex_path.as_posix()
    print(
        f"{rel_path} is missing; restore the tracked algorithm registry source and retry.",
        file=sys.stderr,
    )
    raise SystemExit(1)


def _is_comment_line(path: Path, stripped: str) -> bool:
    if path.suffix == ".py":
        return stripped.startswith("#")
    return stripped.startswith("//")


def parse_code_refs(scan_roots: Sequence[Path] = SCAN_ROOTS) -> dict[str, list[dict]]:
    refs: dict[str, list[dict]] = {}
    for root in scan_roots:
        if not root.exists():
            continue
        for path in sorted(item for item in root.rglob("*") if item.is_file()):
            if path.suffix not in SOURCE_SUFFIXES:
                continue
            lines = path.read_text(encoding="utf-8").splitlines()
            pending: list[tuple[str, int]] = []
            for idx, line in enumerate(lines, start=1):
                stripped = line.strip()
                algid_match = CODE_ALGID_RE.match(stripped)
                if algid_match:
                    pending.append((algid_match.group(1), idx))
                    continue
                if not pending:
                    continue
                if not stripped or _is_comment_line(path, stripped):
                    continue
                try:
                    rel_path = path.relative_to(REPO_ROOT).as_posix()
                except ValueError:
                    rel_path = path.as_posix()
                for algid, comment_line in pending:
                    refs.setdefault(algid, []).append(
                        {
                            "file": rel_path,
                            "line": idx,
                            "comment_line": comment_line,
                            "context": stripped,
                        }
                    )
                pending = []
    return refs


def validate_links(entries: list[dict], code_refs: dict[str, list[dict]]) -> None:
    known_algids = {entry["algid"] for entry in entries}
    unknown = sorted(algid for algid in code_refs if algid not in known_algids)
    if unknown:
        raise ValueError(f"Code AlgID comments reference unknown algorithms: {', '.join(unknown)}")


def attach_code_refs(entries: list[dict], code_refs: dict[str, list[dict]]) -> None:
    for entry in entries:
        entry["owner_refs"] = code_refs.get(entry["algid"], [])


def _status_token(entry: dict) -> str:
    return str(entry.get("status", "")).strip().lower()


def is_documentation_only(entry: dict) -> bool:
    return _status_token(entry) in DOCUMENTATION_ONLY_STATUSES


def missing_owner_entries(entries: list[dict]) -> list[dict]:
    return [entry for entry in entries if not entry.get("owner_refs") and not is_documentation_only(entry)]


def docs_only_entries(entries: list[dict]) -> list[dict]:
    return [entry for entry in entries if is_documentation_only(entry)]


def render_traceability_report(entries: list[dict]) -> str:
    if not entries:
        return "Algorithm traceability: all implementation algorithms have owner comments."

    lines = [
        f"Algorithm traceability warning: {len(entries)} AlgIDs without owner comments:",
    ]
    grouped: dict[str, list[dict]] = defaultdict(list)
    for entry in entries:
        grouped[str(entry.get("section", "<unknown>") or "<unknown>")].append(entry)
    for section in sorted(grouped):
        section_entries = grouped[section]
        lines.append(f"- {section}: {len(section_entries)}")
        for entry in section_entries:
            lines.append(
                f"  - {entry.get('algid', '<unknown>')} "
                f"({entry.get('tex_file', '<unknown>')}:{entry.get('tex_line', '<unknown>')})"
            )
    lines.append("Use a documentation-only or planned status only for explicit non-implementation entries.")
    return "\n".join(lines)


def render_docs_only_audit(entries: list[dict]) -> str:
    docs_entries = docs_only_entries(entries)
    if not docs_entries:
        return "Documentation-only AlgID audit: no documentation-only or planned AlgIDs found."

    lines = [
        f"Documentation-only AlgID audit: {len(docs_entries)} AlgIDs are exempt from strict owner enforcement.",
        "Reason: these entries are documented context or planned routes without a current implementation owner.",
    ]
    grouped: dict[str, list[dict]] = defaultdict(list)
    for entry in docs_entries:
        grouped[str(entry.get("section", "<unknown>") or "<unknown>")].append(entry)
    for section in sorted(grouped):
        section_entries = grouped[section]
        lines.append(f"- {section}: {len(section_entries)}")
        for entry in section_entries:
            description = str(entry.get("description", "") or "").strip()
            suffix = f" - {description}" if description else ""
            lines.append(
                f"  - {entry.get('algid', '<unknown>')} "
                f"({entry.get('tex_file', '<unknown>')}:{entry.get('tex_line', '<unknown>')}){suffix}"
            )
    return "\n".join(lines)


def enforce_traceability(entries: list[dict]) -> None:
    missing_entries = missing_owner_entries(entries)
    if missing_entries:
        raise SystemExit(render_traceability_report(missing_entries))


def yaml_quote(value: object) -> str:
    return json.dumps("" if value is None else value, ensure_ascii=False)


def render_yaml(entries: list[dict]) -> str:
    out: list[str] = []
    out.append("# Generated from docs/latex/algorithms.tex by scripts/docs/sync_algorithm_registry.py")
    out.append("")
    metadata_fields = (
        "section",
        "subsection",
        "subsubsection",
        "family",
        "status",
        "public_api",
        "backend",
        "dependency",
        "derivative_backend",
        "solver_role",
        "implementation_owner",
        "validation",
        "capability_key",
        "description",
        "change_note",
        "tex_file",
    )
    for entry in entries:
        out.append(f"- algid: {yaml_quote(entry['algid'])}")
        for key in metadata_fields:
            out.append(f"  {key}: {yaml_quote(entry.get(key, ''))}")
        out.append(f"  tex_line: {entry['tex_line']}")
        out.append("  body: |")
        body = str(entry.get("body", "")).splitlines() or [""]
        for line in body:
            out.append(f"    {line}" if line else "")
        if entry["owner_refs"]:
            out.append("  owner_refs:")
            for ref in entry["owner_refs"]:
                out.append(f"    - file: {yaml_quote(ref['file'])}")
                out.append(f"      line: {ref['line']}")
                out.append(f"      comment_line: {ref['comment_line']}")
                out.append(f"      context: {yaml_quote(ref['context'])}")
        else:
            out.append("  owner_refs: []")
    out.append("")
    return "\n".join(out).rstrip() + "\n"


def render_latex_body_markdown(body: str) -> str:
    rendered = DISPLAY_MATH_RE.sub(
        lambda match: "\n\n```tex\n" + match.group(1).strip() + "\n```\n\n",
        body,
    )
    rendered = INLINE_MATH_RE.sub(lambda match: f"`{match.group(1).strip()}`", rendered)
    rendered = TEXTTT_RE.sub(lambda match: f"`{match.group(1)}`", rendered)
    return rendered


def body_math_fragments(body: str) -> list[str]:
    fragments: list[str] = []
    fragments.extend(match.group(1).strip() for match in DISPLAY_MATH_RE.finditer(body))
    fragments.extend(match.group(1).strip() for match in INLINE_MATH_RE.finditer(body))
    return [fragment for fragment in fragments if fragment]


def render_markdown(entries: list[dict], rendered_math_map: dict[str, list[str]]) -> str:
    out: list[str] = []
    out.append("# Algorithm Index")
    out.append("")
    out.append("This file is generated from `docs/latex/algorithms.tex` by `scripts/docs/sync_algorithm_registry.py`.")
    out.append(
        "The LaTeX document remains the current source of truth; this Markdown view and `docs/algorithms_registry.yaml` stay aligned with it."
    )
    out.append("")

    current_section = None
    current_subsection = None
    current_subsubsection = None
    for entry in entries:
        section = entry["section"]
        if section != current_section:
            current_section = section
            current_subsection = None
            current_subsubsection = None
            out.append(f"## {section}")
            out.append("")

        subsection = str(entry.get("subsection", "") or "")
        if subsection and subsection != current_subsection:
            current_subsection = subsection
            current_subsubsection = None
            out.append(f"### {subsection}")
            out.append("")

        subsubsection = str(entry.get("subsubsection", "") or "")
        if subsubsection and subsubsection != current_subsubsection:
            current_subsubsection = subsubsection
            out.append(f"#### {subsubsection}")
            out.append("")

        heading_level = 3
        if subsection:
            heading_level += 1
        if subsubsection:
            heading_level += 1
        out.append(f"{'#' * heading_level} `{entry['algid']}`")
        for label, key in (
            ("Family", "family"),
            ("Status", "status"),
            ("Public API", "public_api"),
            ("Backend", "backend"),
            ("Dependency", "dependency"),
            ("Derivative backend", "derivative_backend"),
            ("Solver role", "solver_role"),
            ("Implementation owner", "implementation_owner"),
            ("Validation", "validation"),
            ("Capability key", "capability_key"),
            ("Description", "description"),
            ("Change note", "change_note"),
        ):
            value = str(entry.get(key, "") or "").strip()
            if value:
                out.append(f"- {label}: {value}")
        out.append(f"- LaTeX: `{entry['tex_file']}:{entry['tex_line']}`")
        if entry["owner_refs"]:
            refs = ", ".join(f"`{ref['file']}:{ref['line']}` ({ref['context']})" for ref in entry["owner_refs"])
            out.append(f"- Code owners: {refs}")
        elif is_documentation_only(entry):
            out.append("- Code owners: Documentation-only or planned entry; no current owner expected.")
        else:
            out.append("- Code owners: No `AlgID` owner comment has been attached yet.")
        body = str(entry.get("body", "") or "").strip()
        if body:
            out.append("")
            out.append(render_latex_body_markdown(body))
            out.append("")
            out.append("**LaTeX source**")
            out.append("")
            out.append("```tex")
            out.extend(body.splitlines() or [""])
            out.append("```")
            out.append("")
            rendered_math = rendered_math_map.get(str(entry["algid"]), [])
            if rendered_math:
                out.append("**Rendered formulae**")
                out.append("")
                for fragment in rendered_math:
                    out.append("$$")
                    out.extend(fragment.splitlines() or [""])
                    out.append("$$")
                    out.append("")
        out.append("")

    return "\n".join(out).rstrip() + "\n"


def build_algorithm_math_map(entries: list[dict]) -> dict[str, list[str]]:
    math_map: dict[str, list[str]] = {}
    for entry in entries:
        algid = str(entry["algid"])
        fragments = body_math_fragments(str(entry.get("body", "") or ""))
        math_map[algid] = fragments
    return math_map


def sync_legacy_render_dir(*, check: bool) -> None:
    if LEGACY_RENDER_DIR.exists():
        if check:
            raise SystemExit(
                f"{LEGACY_RENDER_DIR} still contains legacy rendered algorithm images."
            )
        shutil.rmtree(LEGACY_RENDER_DIR)
    if LEGACY_RENDER_ROOT.exists() and not any(LEGACY_RENDER_ROOT.iterdir()):
        LEGACY_RENDER_ROOT.rmdir()


def write_if_changed(path: Path, content: str) -> bool:
    current = path.read_text(encoding="utf-8") if path.exists() else None
    if current == content:
        return False
    path.write_text(content, encoding="utf-8")
    return True


def check_matches(path: Path, expected: str) -> None:
    current = path.read_text(encoding="utf-8") if path.exists() else ""
    if current != expected:
        diff = "".join(
            difflib.unified_diff(
                current.splitlines(keepends=True),
                expected.splitlines(keepends=True),
                fromfile=str(path),
                tofile=f"{path} (expected)",
            )
        )
        raise SystemExit(f"{path} is out of date.\n{diff}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate algorithm registry and Markdown index from docs/latex/algorithms.tex."
    )
    parser.add_argument("--check", action="store_true", help="Validate that generated outputs are up to date.")
    parser.add_argument(
        "--strict-traceability",
        action="store_true",
        help="Fail when any non-documentation-only AlgID lacks an owner comment.",
    )
    parser.add_argument(
        "--docs-only-audit",
        action="store_true",
        help="Print a non-failing audit of documentation-only or planned AlgIDs grouped by section.",
    )
    args = parser.parse_args()

    require_algorithms_tex()
    entries = parse_algorithms(TEX_PATH)
    code_refs = parse_code_refs(SCAN_ROOTS)
    validate_links(entries, code_refs)
    attach_code_refs(entries, code_refs)
    traceability_report = render_traceability_report(missing_owner_entries(entries))
    math_map = build_algorithm_math_map(entries)

    yaml_text = render_yaml(entries)
    markdown_text = render_markdown(entries, math_map)

    if args.check:
        check_matches(REGISTRY_PATH, yaml_text)
        check_matches(MARKDOWN_PATH, markdown_text)
        sync_legacy_render_dir(check=True)
        print("Algorithm registry outputs are up to date.")
        print(traceability_report)
        if args.docs_only_audit:
            print(render_docs_only_audit(entries))
        if args.strict_traceability:
            enforce_traceability(entries)
        return

    if args.strict_traceability:
        enforce_traceability(entries)
    if args.docs_only_audit:
        print(render_docs_only_audit(entries))

    sync_legacy_render_dir(check=False)
    changed_yaml = write_if_changed(REGISTRY_PATH, yaml_text)
    changed_md = write_if_changed(MARKDOWN_PATH, markdown_text)
    if changed_yaml or changed_md:
        print("Updated algorithm registry outputs.")
    else:
        print("Algorithm registry outputs already up to date.")
    print(traceability_report)


if __name__ == "__main__":
    main()
