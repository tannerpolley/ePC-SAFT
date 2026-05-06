from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def asset_roots(repo_root: Path = REPO_ROOT) -> tuple[Path, ...]:
    return (repo_root / "scripts", repo_root / "tests", repo_root / "src")


def should_include_png(path: Path) -> bool:
    parts_lower = tuple(part.lower() for part in path.parts)
    if "__pycache__" in parts_lower:
        return False
    if any(part.startswith("_tmp") for part in parts_lower):
        return False
    if "out" not in parts_lower:
        return False
    return True


def sort_key(path: Path) -> tuple[str, ...]:
    return tuple(part.lower() for part in path.parts)


def collect_pngs(repo_root: Path = REPO_ROOT) -> list[Path]:
    pngs: list[Path] = []
    for asset_root in asset_roots(repo_root):
        if not asset_root.exists():
            continue
        pngs.extend(path for path in asset_root.rglob("out/**/*.png") if should_include_png(path))
    return sorted(
        pngs,
        key=lambda path: sort_key(path.relative_to(repo_root)),
    )


def direct_pngs(root: Path) -> list[Path]:
    if not root.exists():
        return []
    return sorted(
        [path for path in root.glob("*.png") if should_include_png(path)],
        key=lambda path: sort_key(path.relative_to(root)),
    )


def source_path_for_output(rel_path: str) -> str:
    parts = Path(rel_path).parts
    if "out" not in parts:
        return rel_path
    out_index = parts.index("out")
    return Path(*parts[:out_index], *parts[out_index + 1 :]).as_posix()


def folder_for(path: str) -> str:
    folder = str(Path(path).parent).replace("\\", "/")
    return "" if folder == "." else folder


def asset_rows(pngs: list[Path], *, repo_root: Path = REPO_ROOT) -> list[dict[str, str]]:
    rows = []
    for png in pngs:
        output_path = png.relative_to(repo_root).as_posix()
        svg = png.with_suffix(".svg")
        data = png.parent / f"{png.stem}_plot_data.csv"
        rows.append(
            {
                "output_path": output_path,
                "source_path": source_path_for_output(output_path),
                "source_folder": folder_for(source_path_for_output(output_path)),
                "output_folder": folder_for(output_path),
                "svg_path": svg.relative_to(repo_root).as_posix() if svg.exists() else "",
                "data_path": data.relative_to(repo_root).as_posix() if data.exists() else "",
                "name": png.name,
                "title": png.stem.replace("_", " ").replace("-", " "),
            }
        )
    return rows
