from __future__ import annotations

import argparse
import io
from pathlib import Path, PurePosixPath
import xml.etree.ElementTree as ET

REPO_ROOT = Path(__file__).resolve().parents[2]
IDEA_DIR = REPO_ROOT / ".idea"
MODULE_URL_PREFIX = "file://$MODULE_DIR$"
CONTENT_URL = MODULE_URL_PREFIX
TRANSIENT_PATHS: tuple[str, ...] = (
    "build",
    "dist",
    ".venv",
    ".worktrees",
    "_codex",
    ".pytest_cache",
    ".ruff_cache",
    ".mypy_cache",
    "docs/_build",
    "docs/latex/out",
    "results/runs",
)
UNIVERSAL_SOURCE_ROOTS: tuple[tuple[str, bool], ...] = (
    ("src", False),
    ("tests", True),
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Normalize JetBrains module metadata under .idea/.",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="Report pending changes without writing files.",
    )
    mode.add_argument(
        "--apply",
        action="store_true",
        help="Rewrite module files in place.",
    )
    return parser


def _module_url(relative_path: str) -> str:
    normalized = relative_path.replace("\\", "/").strip("/")
    if not normalized:
        return CONTENT_URL
    return f"{MODULE_URL_PREFIX}/{normalized}"


def _relative_module_path(url: str | None) -> PurePosixPath | None:
    if url is None:
        return None
    if url == MODULE_URL_PREFIX:
        return PurePosixPath(".")
    prefix = f"{MODULE_URL_PREFIX}/"
    if not url.startswith(prefix):
        return None
    relative = url[len(prefix) :].strip("/")
    if not relative:
        return PurePosixPath(".")
    return PurePosixPath(relative)


def _is_under_path(url: str | None, relative_path: str) -> bool:
    candidate = _relative_module_path(url)
    if candidate is None:
        return False
    base = PurePosixPath(relative_path)
    return candidate == base or base in candidate.parents


def _existing_transient_paths() -> tuple[str, ...]:
    existing: list[str] = []
    for relative_path in TRANSIENT_PATHS:
        if (REPO_ROOT / relative_path).exists():
            existing.append(relative_path)
    return tuple(existing)


def _discover_iml_files() -> tuple[Path, ...]:
    candidates = list(IDEA_DIR.glob("*.iml"))
    candidates.extend(REPO_ROOT.glob("*.iml"))
    return tuple(sorted(path.resolve() for path in candidates))


def _find_child(parent: ET.Element, tag: str, **attrs: str) -> ET.Element | None:
    for child in parent.findall(tag):
        if all(child.get(key) == value for key, value in attrs.items()):
            return child
    return None


def _ensure_root_manager(module_root: ET.Element, actions: list[str]) -> ET.Element:
    manager = _find_child(module_root, "component", name="NewModuleRootManager")
    if manager is None:
        manager = ET.Element("component", {"name": "NewModuleRootManager", "inherit-compiler-output": "true"})
        module_root.insert(0, manager)
        actions.append("create NewModuleRootManager component")
    elif manager.get("inherit-compiler-output") is None:
        manager.set("inherit-compiler-output", "true")
        actions.append("set NewModuleRootManager inherit-compiler-output=true")
    return manager


def _ensure_exclude_output(manager: ET.Element, actions: list[str]) -> None:
    exclude_output = manager.find("exclude-output")
    if exclude_output is None:
        manager.insert(0, ET.Element("exclude-output"))
        actions.append("add exclude-output element")


def _ensure_content_root(manager: ET.Element, actions: list[str]) -> ET.Element:
    content = _find_child(manager, "content", url=CONTENT_URL)
    if content is not None:
        return content

    first_content = manager.find("content")
    if first_content is not None:
        if first_content.get("url") != CONTENT_URL:
            first_content.set("url", CONTENT_URL)
            actions.append(f"set content root to {CONTENT_URL}")
        return first_content

    content = ET.Element("content", {"url": CONTENT_URL})
    insert_at = 1 if manager.find("exclude-output") is not None else 0
    manager.insert(insert_at, content)
    actions.append(f"create content root {CONTENT_URL}")
    return content


def _ensure_source_folder_order_entry(manager: ET.Element, actions: list[str]) -> None:
    for entry in manager.findall("orderEntry"):
        if entry.get("type") == "sourceFolder":
            return
    manager.append(ET.Element("orderEntry", {"type": "sourceFolder", "forTests": "false"}))
    actions.append("add sourceFolder orderEntry")


def _prune_stale_module_dependencies(manager: ET.Element, declared_modules: set[str], actions: list[str]) -> None:
    if not declared_modules:
        return
    for entry in list(manager.findall("orderEntry")):
        if entry.get("type") != "module":
            continue
        module_name = entry.get("module-name")
        if module_name and module_name not in declared_modules:
            manager.remove(entry)
            actions.append(f"remove stale module dependency {module_name}")


def _replace_content_roots(content: ET.Element, transient_paths: tuple[str, ...], actions: list[str]) -> None:
    existing_source_folders = list(content.findall("sourceFolder"))
    existing_exclude_folders = list(content.findall("excludeFolder"))
    preserved_children = [child for child in list(content) if child.tag not in {"sourceFolder", "excludeFolder"}]

    desired_sources: dict[tuple[str, str], ET.Element] = {}
    for source in existing_source_folders:
        url = source.get("url")
        is_test = "true" if source.get("isTestSource") == "true" else "false"
        if any(_is_under_path(url, path) for path in transient_paths):
            actions.append(f"remove transient sourceFolder {url}")
            continue
        key = (url or "", is_test)
        if key in desired_sources:
            actions.append(f"remove duplicate sourceFolder {url} (isTestSource={is_test})")
            continue
        desired_sources[key] = ET.Element(
            "sourceFolder",
            {"url": url or "", "isTestSource": is_test},
        )

    for relative_path, is_test in UNIVERSAL_SOURCE_ROOTS:
        url = _module_url(relative_path)
        key = (url, "true" if is_test else "false")
        if key not in desired_sources:
            desired_sources[key] = ET.Element(
                "sourceFolder",
                {"url": url, "isTestSource": "true" if is_test else "false"},
            )
            label = "test source root" if is_test else "source root"
            actions.append(f"ensure {label} {url}")

    desired_excludes: dict[str, ET.Element] = {}
    for exclude in existing_exclude_folders:
        url = exclude.get("url")
        if not url:
            continue
        nested_transient_parent = next(
            (
                relative_path
                for relative_path in TRANSIENT_PATHS
                if _is_under_path(url, relative_path) and url != _module_url(relative_path)
            ),
            None,
        )
        if nested_transient_parent is not None:
            actions.append(f"remove nested transient excludeFolder {url}")
            continue
        if url not in desired_excludes:
            desired_excludes[url] = ET.Element("excludeFolder", {"url": url})
        else:
            actions.append(f"remove duplicate excludeFolder {url}")

    for relative_path in transient_paths:
        url = _module_url(relative_path)
        if url not in desired_excludes:
            desired_excludes[url] = ET.Element("excludeFolder", {"url": url})
            actions.append(f"ensure excludeFolder {url}")

    sorted_exclude_urls = sorted(desired_excludes)
    collapsed_excludes: list[ET.Element] = []
    kept_exclude_urls: list[str] = []
    for url in sorted_exclude_urls:
        if any(_is_under_path(url, kept_url[len(MODULE_URL_PREFIX) + 1 :]) for kept_url in kept_exclude_urls if kept_url.startswith(f"{MODULE_URL_PREFIX}/")):
            actions.append(f"remove redundant nested excludeFolder {url}")
            continue
        collapsed_excludes.append(desired_excludes[url])
        kept_exclude_urls.append(url)

    sorted_sources = sorted(
        desired_sources.values(),
        key=lambda element: (
            element.get("url", ""),
            element.get("isTestSource", "false") == "true",
        ),
    )

    for child in list(content):
        content.remove(child)
    for source in sorted_sources:
        content.append(source)
    for exclude in collapsed_excludes:
        content.append(exclude)
    for child in preserved_children:
        content.append(child)


def _serialize_tree(tree: ET.ElementTree) -> str:
    ET.indent(tree, space="  ")
    buffer = io.BytesIO()
    tree.write(buffer, encoding="UTF-8", xml_declaration=True)
    return buffer.getvalue().decode("UTF-8")


def _load_declared_modules() -> set[str]:
    modules_path = IDEA_DIR / "modules.xml"
    if not modules_path.exists():
        return set()
    tree = ET.parse(modules_path)
    root = tree.getroot()
    names: set[str] = set()
    for module in root.findall(".//module"):
        for attr_name in ("filepath", "fileurl"):
            value = module.get(attr_name)
            if not value or not value.endswith(".iml"):
                continue
            names.add(Path(value).stem)
            break
    return names


def _module_dependency_warnings(iml_path: Path, module_root: ET.Element, declared_modules: set[str]) -> list[str]:
    if not declared_modules:
        return []
    warnings: list[str] = []
    for dependency in module_root.findall(".//orderEntry[@type='module']"):
        module_name = dependency.get("module-name")
        if module_name and module_name not in declared_modules:
            warnings.append(
                f"stale module dependency '{module_name}' is not declared in .idea/modules.xml"
            )
    return warnings


def _normalize_iml(iml_path: Path, transient_paths: tuple[str, ...], declared_modules: set[str]) -> tuple[list[str], list[str], str | None]:
    original_text = iml_path.read_text(encoding="UTF-8")
    tree = ET.parse(iml_path)
    module_root = tree.getroot()
    actions: list[str] = []
    manager = _ensure_root_manager(module_root, actions)
    _ensure_exclude_output(manager, actions)
    content = _ensure_content_root(manager, actions)
    _ensure_source_folder_order_entry(manager, actions)
    _prune_stale_module_dependencies(manager, declared_modules, actions)
    _replace_content_roots(content, transient_paths, actions)
    warnings = _module_dependency_warnings(iml_path, module_root, declared_modules)
    proposed_text = _serialize_tree(tree)
    if proposed_text == original_text:
        actions = []
    return actions, warnings, proposed_text if proposed_text != original_text else None


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    transient_paths = _existing_transient_paths()
    declared_modules = _load_declared_modules()
    iml_files = _discover_iml_files()

    pending_changes = 0
    warnings_found = 0
    for iml_path in iml_files:
        actions, warnings, proposed_text = _normalize_iml(iml_path, transient_paths, declared_modules)
        relative_path = iml_path.relative_to(REPO_ROOT).as_posix()
        for warning in warnings:
            warnings_found += 1
            print(f"WARNING {relative_path}: {warning}")
        if not actions:
            print(f"OK {relative_path}: no changes")
            continue

        pending_changes += 1
        prefix = "APPLY" if args.apply else "DRY-RUN"
        for action in actions:
            print(f"{prefix} {relative_path}: {action}")
        if args.apply and proposed_text is not None:
            iml_path.write_text(proposed_text, encoding="UTF-8", newline="\n")

    if pending_changes == 0:
        print("No pending module metadata changes.")
    if warnings_found == 0:
        print("No stale module dependency warnings.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
