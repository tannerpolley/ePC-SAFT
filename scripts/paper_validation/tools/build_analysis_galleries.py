from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
PLOTS_ROOT = REPO_ROOT / "docs" / "plots"
MANIFEST_PATH = PLOTS_ROOT / "manifest.json"


def gallery_repo_root() -> Path:
    if PLOTS_ROOT.name == "plots" and PLOTS_ROOT.parent.name == "docs":
        return PLOTS_ROOT.parent.parent
    return REPO_ROOT


def asset_roots() -> tuple[Path, ...]:
    root = gallery_repo_root()
    return (root / "scripts", root / "tests", root / "src")


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


def collect_pngs(root: Path) -> list[Path]:
    roots = asset_roots() if root == PLOTS_ROOT else (root,)
    pngs: list[Path] = []
    for asset_root in roots:
        if not asset_root.exists():
            continue
        pngs.extend(path for path in asset_root.rglob("out/**/*.png") if should_include_png(path))
    return sorted(
        pngs,
        key=lambda path: sort_key(path.relative_to(REPO_ROOT)),
    )


def direct_pngs(root: Path) -> list[Path]:
    if not root.exists():
        return []
    return sorted(
        [path for path in root.glob("*.png") if should_include_png(path)],
        key=lambda path: sort_key(path.relative_to(root)),
    )


def _source_path_for_output(rel_path: str) -> str:
    parts = Path(rel_path).parts
    if "out" not in parts:
        return rel_path
    out_index = parts.index("out")
    return Path(*parts[:out_index], *parts[out_index + 1 :]).as_posix()


def _folder_for(path: str) -> str:
    folder = str(Path(path).parent).replace("\\", "/")
    return "" if folder == "." else folder


def image_manifest(pngs: list[Path]) -> list[dict[str, str]]:
    manifest = []
    for png in pngs:
        repo_root = gallery_repo_root()
        output_path = png.relative_to(repo_root).as_posix()
        asset_path = Path("../../") / output_path
        svg = png.with_suffix(".svg")
        data = png.parent / f"{png.stem}_plot_data.csv"
        source_path = _source_path_for_output(output_path)
        manifest.append(
            {
                "path": asset_path.as_posix(),
                "folder": _folder_for(source_path),
                "output_path": output_path,
                "output_folder": _folder_for(output_path),
                "svg_path": (Path("../../") / svg.relative_to(repo_root)).as_posix() if svg.exists() else "",
                "data_path": (Path("../../") / data.relative_to(repo_root)).as_posix() if data.exists() else "",
                "source_path": source_path,
                "source_folder": _folder_for(source_path),
                "name": png.name,
                "title": png.stem.replace("_", " ").replace("-", " "),
            }
        )
    return manifest


_image_manifest = image_manifest


def manifest_path(root: Path = PLOTS_ROOT) -> Path:
    return root / "manifest.json"


def load_manifest(path: Path | None = None, root: Path = PLOTS_ROOT) -> list[dict[str, str]]:
    path = path or manifest_path(root)
    if not path.exists():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    items = payload.get("items", [])
    if not isinstance(items, list):
        return []
    return [dict(item) for item in items if isinstance(item, dict)]


def gallery_entries(root: Path = PLOTS_ROOT, pngs: list[Path] | None = None) -> list[dict[str, str]]:
    if pngs is not None:
        return image_manifest(pngs)
    entries = load_manifest(root=root)
    if entries:
        return entries
    return image_manifest(collect_pngs(root))
