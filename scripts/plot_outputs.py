from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PLOTS_ROOT = REPO_ROOT / "docs" / "plots"
PAPER_VALIDATION_SOURCE_ROOT = REPO_ROOT / "scripts" / "paper_validation"
PAPER_VALIDATION_PLOTS_ROOT = PLOTS_ROOT / "paper_validation"
FITS_PLOTS_ROOT = PLOTS_ROOT / "fits"


def _clean_analysis_name(name: str) -> str:
    return name.removesuffix("_analysis")


def paper_validation_path(source_path: str | Path, filename: str | None = None) -> Path:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    rel_dir = source_dir.relative_to(PAPER_VALIDATION_SOURCE_ROOT)
    rel_parts = list(rel_dir.parts)
    if rel_parts:
        rel_parts[0] = _clean_analysis_name(rel_parts[0])
    target_dir = PAPER_VALIDATION_PLOTS_ROOT.joinpath(*rel_parts)
    target = target_dir / (filename if filename is not None else source.name)
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def paper_validation_dir(source_path: str | Path) -> Path:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    rel_dir = source_dir.relative_to(PAPER_VALIDATION_SOURCE_ROOT)
    rel_parts = list(rel_dir.parts)
    if rel_parts:
        rel_parts[0] = _clean_analysis_name(rel_parts[0])
    target = PAPER_VALIDATION_PLOTS_ROOT.joinpath(*rel_parts)
    target.mkdir(parents=True, exist_ok=True)
    return target


def paper_validation_output_path(path: str | Path) -> Path:
    source = Path(path).resolve()
    if source.is_relative_to(PAPER_VALIDATION_SOURCE_ROOT):
        return paper_validation_path(source.parent, source.name)
    source.parent.mkdir(parents=True, exist_ok=True)
    return source


def fits_plot_path(*parts: str | Path) -> Path:
    target = FITS_PLOTS_ROOT.joinpath(*(str(part) for part in parts))
    target.parent.mkdir(parents=True, exist_ok=True)
    return target
