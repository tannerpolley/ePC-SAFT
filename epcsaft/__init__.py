from __future__ import annotations

from pathlib import Path

_SRC_PACKAGE = Path(__file__).resolve().parents[1] / "src" / "epcsaft"
if not _SRC_PACKAGE.is_dir():
    raise ModuleNotFoundError("Cannot find source package at src/epcsaft.")

__path__ = [str(_SRC_PACKAGE)]
__file__ = str(_SRC_PACKAGE / "__init__.py")

with open(__file__, "rb") as _handle:
    exec(compile(_handle.read(), __file__, "exec"), globals())
