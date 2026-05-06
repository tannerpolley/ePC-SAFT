from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.paper_validation.tools import build_analysis_galleries

DEFAULT_MANIFEST = REPO_ROOT / "docs" / "plots" / "manifest.json"
MANIFEST_VERSION = 1
REQUIRED_FIELDS = {
    "path",
    "folder",
    "output_path",
    "output_folder",
    "svg_path",
    "data_path",
    "source_path",
    "source_folder",
    "name",
    "title",
}

CANONICAL_TEXT_REPLACEMENTS = (
    ("maic_m_", "miac_m_"),
    ("maic m", "miac m"),
    ("diagnostics/output/plots", "diagnostics/out/plots"),
    ("diagnostics/output/plots/", "diagnostics/out/plots/"),
    ("diagnostics/out/plots/out/", "diagnostics/out/plots/"),
    ("diagnostics/out/plots/figure6b_fit_checks/out/", "diagnostics/out/plots/figure6b_fit_checks/"),
    (
        "diagnostics/out/plots/figure6b_fit_checks_best_scalar_zfits/out/",
        "diagnostics/out/plots/figure6b_fit_checks_best_scalar_zfits/",
    ),
)


def canonicalize_item(item: dict[str, str]) -> dict[str, str]:
    canonical = dict(item)
    for key, value in canonical.items():
        if not isinstance(value, str):
            continue
        for old, new in CANONICAL_TEXT_REPLACEMENTS:
            value = value.replace(old, new)
        canonical[key] = value
    return canonical


def canonicalize_items(items: list[dict[str, str]]) -> list[dict[str, str]]:
    return [canonicalize_item(item) for item in items]


def manifest_payload(items: list[dict[str, str]]) -> dict[str, Any]:
    return {
        "version": MANIFEST_VERSION,
        "description": "Static plot metadata for source-owned ePC-SAFT plot outputs.",
        "items": canonicalize_items(items),
    }


def collect_manifest_items() -> list[dict[str, str]]:
    pngs = build_analysis_galleries.collect_pngs(build_analysis_galleries.PLOTS_ROOT)
    return build_analysis_galleries.image_manifest(pngs)


def read_manifest(path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_manifest(path: Path = DEFAULT_MANIFEST) -> Path:
    payload = manifest_payload(collect_manifest_items())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def validate_manifest(path: Path = DEFAULT_MANIFEST) -> list[str]:
    errors: list[str] = []
    if not path.exists():
        return [f"Missing plot manifest: {path}"]
    try:
        payload = read_manifest(path)
    except json.JSONDecodeError as exc:
        return [f"Invalid JSON in {path}: {exc}"]

    if payload.get("version") != MANIFEST_VERSION:
        errors.append(f"Manifest version must be {MANIFEST_VERSION}.")
    items = payload.get("items")
    if not isinstance(items, list):
        errors.append("Manifest must contain an items list.")
        return errors

    seen: set[str] = set()
    for index, item in enumerate(items):
        if not isinstance(item, dict):
            errors.append(f"Item {index} is not an object.")
            continue
        missing = sorted(REQUIRED_FIELDS - set(item))
        if missing:
            errors.append(f"Item {index} is missing fields: {', '.join(missing)}")
        output_path = str(item.get("output_path", ""))
        if not output_path:
            errors.append(f"Item {index} has an empty output_path.")
        if output_path in seen:
            errors.append(f"Duplicate output_path: {output_path}")
        seen.add(output_path)
        if ".html" in output_path.lower():
            errors.append(f"Manifest item should not reference HTML assets: {output_path}")
    return errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build or validate the tracked plot manifest.")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="Manifest JSON path.")
    parser.add_argument(
        "--refresh", action="store_true", help="Refresh manifest from currently generated local PNG assets."
    )
    parser.add_argument(
        "--check", action="store_true", help="Validate the tracked manifest without requiring generated assets."
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.refresh:
        path = write_manifest(args.manifest)
        count = len(read_manifest(path).get("items", []))
        print(f"Wrote {path} with {count} item(s).")
        return 0

    errors = validate_manifest(args.manifest)
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print(f"Plot manifest is valid: {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
