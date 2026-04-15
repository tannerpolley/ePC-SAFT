from __future__ import annotations

import argparse
import difflib
import json
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TEX_PATH = REPO_ROOT / "docs" / "latex" / "equations.tex"
MARKDOWN_PATH = REPO_ROOT / "docs" / "equations.md"
REGISTRY_PATH = REPO_ROOT / "docs" / "equations_registry.yaml"
NATIVE_ROOT = REPO_ROOT / "src" / "epcsaft" / "native"

SECTION_RE = re.compile(r"\\section\{(.+?)\}")
EQID_RE = re.compile(r"%\s*EqID:\s*([A-Za-z0-9_]+)")
META_RE = re.compile(r"%\s*([A-Za-z][A-Za-z ]+):\s*(.*)")
LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
CODE_EQID_RE = re.compile(r"//\s*EqID:\s*([A-Za-z0-9_]+)")
BEGIN_ENV_RE = re.compile(r"\\begin\{([A-Za-z*]+)\}")


def parse_equations(tex_path: Path) -> list[dict]:
    lines = tex_path.read_text(encoding="utf-8").splitlines()
    entries: list[dict] = []
    section = ""
    i = 0
    seen_eqids: set[str] = set()
    while i < len(lines):
        line = lines[i]
        section_match = SECTION_RE.search(line)
        if section_match:
            section = section_match.group(1).strip()
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
    out.append("The LaTeX document remains the current source of truth; this Markdown view and `docs/equations_registry.yaml` stay aligned with it.")
    out.append("")

    current_section = None
    for entry in entries:
        section = entry["section"]
        if section != current_section:
            current_section = section
            out.append(f"## {section}")
            out.append("")

        out.append(f"### `{entry['eqid']}`")
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
            refs = ", ".join(
                f"`{ref['file']}:{ref['line']}` ({ref['context']})"
                for ref in entry["cpp_refs"]
            )
            out.append(f"- C++: {refs}")
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
    parser = argparse.ArgumentParser(description="Generate equation registry and Markdown index from docs/latex/equations.tex.")
    parser.add_argument("--check", action="store_true", help="Validate that generated outputs are up to date.")
    args = parser.parse_args()

    entries = parse_equations(TEX_PATH)
    code_refs = parse_code_refs(NATIVE_ROOT)
    validate_links(entries, code_refs)
    attach_code_refs(entries, code_refs)

    yaml_text = render_yaml(entries)
    markdown_text = render_markdown(entries)

    if args.check:
        check_matches(REGISTRY_PATH, yaml_text)
        check_matches(MARKDOWN_PATH, markdown_text)
        print("Equation registry outputs are up to date.")
        return

    changed_yaml = write_if_changed(REGISTRY_PATH, yaml_text)
    changed_md = write_if_changed(MARKDOWN_PATH, markdown_text)
    if changed_yaml or changed_md:
        print("Updated equation registry outputs.")
    else:
        print("Equation registry outputs already up to date.")


if __name__ == "__main__":
    main()
