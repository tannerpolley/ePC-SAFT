from __future__ import annotations


def require_epcsaft_install() -> None:
    try:
        import epcsaft  # noqa: F401
        from epcsaft import ePCSAFTMixture  # noqa: F401
        from epcsaft import ePCSAFTState  # noqa: F401
        import epcsaft.parameters  # noqa: F401
    except Exception as exc:
        raise RuntimeError(
            "epcsaft must be importable from the active environment with the OOP API. "
            "Run `python scripts/build_epcsaft.py` in the ePC-SAFT environment, then retry."
        ) from exc


