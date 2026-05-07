from __future__ import annotations

import argparse
from collections import defaultdict
import difflib
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
TEX_PATH = REPO_ROOT / "docs" / "latex" / "equations.tex"
MARKDOWN_PATH = REPO_ROOT / "docs" / "equations.md"
REGISTRY_PATH = REPO_ROOT / "docs" / "equations_registry.yaml"
NATIVE_ROOT = REPO_ROOT / "src" / "epcsaft" / "native"

SECTION_RE = re.compile(r"\\section\{(.+?)\}")
SUBSECTION_RE = re.compile(r"\\subsection\{(.+?)\}")
SUBSUBSECTION_RE = re.compile(r"\\subsubsection\{(.+?)\}")
EQID_RE = re.compile(r"%\s*EqID:\s*([A-Za-z0-9_]+)")
META_RE = re.compile(r"%\s*([A-Za-z][A-Za-z ]+):\s*(.*)")
LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
CODE_EQID_RE = re.compile(r"//\s*EqID:\s*([A-Za-z0-9_]+)")
BEGIN_ENV_RE = re.compile(r"\\begin\{([A-Za-z*]+)\}")
DOCUMENTATION_ONLY_STATUSES = {
    "documentation-only",
    "documentation_only",
    "docs-only",
    "docs_only",
    "reference-only",
    "reference_only",
}


def parse_equations(tex_path: Path) -> list[dict]:
    lines = tex_path.read_text(encoding="utf-8").splitlines()
    entries: list[dict] = []
    section = ""
    subsection = ""
    subsubsection = ""
    i = 0
    seen_eqids: set[str] = set()
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

        eqid_match = EQID_RE.match(line.strip())
        if not eqid_match:
            i += 1
            continue

        eqid = eqid_match.group(1)
        if eqid in seen_eqids:
            raise ValueError(f"Duplicate EqID in {tex_path}: {eqid}")
        seen_eqids.add(eqid)

        entry: dict[str, object] = {
            "eqid": eqid,
            "section": section,
            "subsection": subsection,
            "subsubsection": subsubsection,
            "tex_file": tex_path.relative_to(REPO_ROOT).as_posix(),
            "tex_line": i + 1,
        }

        i += 1
        while i < len(lines):
            stripped = lines[i].strip()
            if stripped == "":
                i += 1
                continue
            if BEGIN_ENV_RE.match(stripped):
                break
            meta_match = META_RE.match(stripped)
            if meta_match:
                key = meta_match.group(1).strip().lower().replace(" ", "_")
                entry[key] = meta_match.group(2).strip()
            i += 1

        begin_match = BEGIN_ENV_RE.match(lines[i].strip()) if i < len(lines) else None
        if begin_match is None:
            raise ValueError(f"EqID {eqid} is not followed by an equation block in {tex_path}")
        env_name = begin_match.group(1)
        end_token = rf"\end{{{env_name}}}"

        equation_lines: list[str] = []
        label = ""
        i += 1
        while i < len(lines):
            stripped = lines[i].strip()
            if stripped == end_token:
                break
            label_match = LABEL_RE.search(stripped)
            if label_match:
                label = label_match.group(1)
            elif stripped:
                equation_lines.append(lines[i].rstrip())
            i += 1

        if i >= len(lines):
            raise ValueError(f"Equation block for EqID {eqid} is not terminated in {tex_path}")

        entry["label"] = label
        entry["latex"] = "\n".join(equation_lines).strip()
        entry["cpp_refs"] = []
        entries.append(entry)
        i += 1

    return entries


def require_equations_tex(tex_path: Path = TEX_PATH) -> None:
    """Exit with an actionable message when the tracked LaTeX source is absent."""
    if tex_path.exists():
        return
    try:
        rel_path = tex_path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        rel_path = tex_path.as_posix()
    print(
        "{} is missing; docs/latex is tracked repo content. Restore the file from git or refresh the checkout and retry.".format(
            rel_path,
        ),
        file=sys.stderr,
    )
    raise SystemExit(1)


def parse_code_refs(native_root: Path) -> dict[str, list[dict]]:
    refs: dict[str, list[dict]] = {}
    for path in sorted(native_root.rglob("*")):
        if path.suffix not in {".cpp", ".h"}:
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        pending: list[tuple[str, int]] = []
        for idx, line in enumerate(lines, start=1):
            stripped = line.strip()
            eqid_match = CODE_EQID_RE.match(stripped)
            if eqid_match:
                pending.append((eqid_match.group(1), idx))
                continue
            if not pending:
                continue
            if not stripped or stripped.startswith("//"):
                continue
            context = stripped
            rel_path = path.relative_to(REPO_ROOT).as_posix()
            for eqid, comment_line in pending:
                refs.setdefault(eqid, []).append(
                    {
                        "file": rel_path,
                        "line": idx,
                        "comment_line": comment_line,
                        "context": context,
                    }
                )
            pending = []
    return refs


def validate_links(entries: list[dict], code_refs: dict[str, list[dict]]) -> None:
    known_eqids = {entry["eqid"] for entry in entries}
    unknown = sorted(eqid for eqid in code_refs if eqid not in known_eqids)
    if unknown:
        raise ValueError(f"C++ EqID comments reference unknown equations: {', '.join(unknown)}")


def attach_code_refs(entries: list[dict], code_refs: dict[str, list[dict]]) -> None:
    for entry in entries:
        entry["cpp_refs"] = code_refs.get(entry["eqid"], [])


def _status_token(entry: dict) -> str:
    return str(entry.get("status", "")).strip().lower()


def is_documentation_only(entry: dict) -> bool:
    return _status_token(entry) in DOCUMENTATION_ONLY_STATUSES


def missing_cpp_ref_entries(entries: list[dict]) -> list[dict]:
    return [entry for entry in entries if not entry.get("cpp_refs") and not is_documentation_only(entry)]


def docs_only_entries(entries: list[dict]) -> list[dict]:
    return [entry for entry in entries if is_documentation_only(entry)]


def render_traceability_report(entries: list[dict]) -> str:
    if not entries:
        return "Equation traceability: all implementation equations have C++ owner comments."

    lines = [
        f"Equation traceability warning: {len(entries)} EqIDs without C++ owner comments:",
    ]
    grouped: dict[str, list[dict]] = defaultdict(list)
    for entry in entries:
        grouped[str(entry.get("section", "<unknown>") or "<unknown>")].append(entry)
    for section in sorted(grouped):
        section_entries = grouped[section]
        lines.append(f"- {section}: {len(section_entries)}")
        for entry in section_entries:
            lines.append(
                f"  - {entry.get('eqid', '<unknown>')} ({entry.get('tex_file', '<unknown>')}:{entry.get('tex_line', '<unknown>')})"
            )
    lines.append("Mark documentation-only equations with status: Documentation-only.")
    return "\n".join(lines)


def render_docs_only_audit(entries: list[dict]) -> str:
    docs_entries = docs_only_entries(entries)
    if not docs_entries:
        return "Documentation-only EqID audit: no documentation-only EqIDs found."

    lines = [
        f"Documentation-only EqID audit: {len(docs_entries)} EqIDs are exempt from strict C++ owner enforcement.",
        "Reason: these equations are reference material, notation, derivation, or explanatory context with no direct native owner.",
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
                f"  - {entry.get('eqid', '<unknown>')} "
                f"({entry.get('tex_file', '<unknown>')}:{entry.get('tex_line', '<unknown>')}){suffix}"
            )
    return "\n".join(lines)


def enforce_traceability(entries: list[dict]) -> None:
    missing_entries = missing_cpp_ref_entries(entries)
    if missing_entries:
        raise SystemExit(render_traceability_report(missing_entries))


def yaml_quote(value: object) -> str:
    return json.dumps("" if value is None else value, ensure_ascii=False)


def render_yaml(entries: list[dict]) -> str:
    out: list[str] = []
    out.append("# Generated from docs/latex/equations.tex by scripts/sync_equation_registry.py")
    out.append("")
    for entry in entries:
        out.append(f"- eqid: {yaml_quote(entry['eqid'])}")
        for key in (
            "section",
            "subsection",
            "subsubsection",
            "label",
            "source",
            "status",
            "description",
            "change_note",
            "tex_file",
        ):
            out.append(f"  {key}: {yaml_quote(entry.get(key, ''))}")
        out.append(f"  tex_line: {entry['tex_line']}")
        out.append("  latex: |")
        latex = str(entry.get("latex", "")).splitlines() or [""]
        for line in latex:
            out.append(f"    {line}")
        if entry["cpp_refs"]:
            out.append("  cpp_refs:")
            for ref in entry["cpp_refs"]:
                out.append(f"    - file: {yaml_quote(ref['file'])}")
                out.append(f"      line: {ref['line']}")
                out.append(f"      comment_line: {ref['comment_line']}")
                out.append(f"      context: {yaml_quote(ref['context'])}")
        else:
            out.append("  cpp_refs: []")
    out.append("")
    return "\n".join(out)


def render_markdown(entries: list[dict]) -> str:
    out: list[str] = []
    out.append("# Equation Index")
    out.append("")
    out.append("This file is generated from `docs/latex/equations.tex` by `scripts/sync_equation_registry.py`.")
    out.append(
        "The LaTeX document remains the current source of truth; this Markdown view and `docs/equations_registry.yaml` stay aligned with it."
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
        out.append(f"{'#' * heading_level} `{entry['eqid']}`")
        out.append(f"- Label: `{entry.get('label', '')}`")
        if entry.get("source"):
            out.append(f"- Source: {entry['source']}")
        if entry.get("status"):
            out.append(f"- Status: {entry['status']}")
        if entry.get("description"):
            out.append(f"- Description: {entry['description']}")
        if entry.get("change_note"):
            out.append(f"- Change note: {entry['change_note']}")
        out.append(f"- LaTeX: `{entry['tex_file']}:{entry['tex_line']}`")
        if entry["cpp_refs"]:
            refs = ", ".join(f"`{ref['file']}:{ref['line']}` ({ref['context']})" for ref in entry["cpp_refs"])
            out.append(f"- C++: {refs}")
        elif is_documentation_only(entry):
            out.append("- C++: Documentation-only: no direct native owner expected.")
        else:
            out.append("- C++: No `EqID` owner comment has been attached yet.")
        out.append("")
        out.append("```tex")
        out.extend(str(entry.get("latex", "")).splitlines() or [""])
        out.append("```")
        out.append("")

    return "\n".join(out)


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
        description="Generate equation registry and Markdown index from docs/latex/equations.tex."
    )
    parser.add_argument("--check", action="store_true", help="Validate that generated outputs are up to date.")
    parser.add_argument(
        "--strict-traceability",
        action="store_true",
        help="Fail when any non-documentation-only EqID lacks a C++ owner comment.",
    )
    parser.add_argument(
        "--docs-only-audit",
        action="store_true",
        help="Print a non-failing audit of documentation-only EqIDs grouped by section.",
    )
    args = parser.parse_args()

    require_equations_tex()
    entries = parse_equations(TEX_PATH)
    code_refs = parse_code_refs(NATIVE_ROOT)
    validate_links(entries, code_refs)
    attach_code_refs(entries, code_refs)
    traceability_report = render_traceability_report(missing_cpp_ref_entries(entries))

    yaml_text = render_yaml(entries)
    markdown_text = render_markdown(entries)

    if args.check:
        check_matches(REGISTRY_PATH, yaml_text)
        check_matches(MARKDOWN_PATH, markdown_text)
        print("Equation registry outputs are up to date.")
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

    changed_yaml = write_if_changed(REGISTRY_PATH, yaml_text)
    changed_md = write_if_changed(MARKDOWN_PATH, markdown_text)
    if changed_yaml or changed_md:
        print("Updated equation registry outputs.")
    else:
        print("Equation registry outputs already up to date.")
    print(traceability_report)


if __name__ == "__main__":
    main()
