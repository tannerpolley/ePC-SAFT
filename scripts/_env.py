from __future__ import annotations


def require_pcsaft_install() -> None:
    try:
        import pcsaft  # noqa: F401
        from pcsaft import PCSAFTMixture  # noqa: F401
        from pcsaft import PCSAFTState  # noqa: F401
        import pcsaft.parameters  # noqa: F401
    except Exception as exc:
        raise RuntimeError(
            "pcsaft must be importable from the active environment with the OOP API. "
            "Run `python scripts/build_pcsaft.py` in the PC-SAFT environment, then retry."
        ) from exc
