from __future__ import annotations


def require_epcsaft_install() -> None:
    try:
        import epcsaft
        import epcsaft.parameters
        from epcsaft import (
            ePCSAFTMixture,
            ePCSAFTState,
        )
    except Exception as exc:
        raise RuntimeError(
            "epcsaft must be importable from the active environment with the OOP API. "
            "Run `uv sync --no-install-project`, then `uv run python scripts/build_epcsaft.py`, then retry."
        ) from exc
