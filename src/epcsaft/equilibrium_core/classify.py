"""Route classification for Python equilibrium calculations."""

from __future__ import annotations

from typing import Any

import numpy as np


def classify_equilibrium_route(mixture: Any, kind: str, backend: str | None = None) -> dict[str, str]:
    """Return the internal route for a public equilibrium request."""
    token = str(kind).strip().lower()
    charges = np.asarray(mixture.parameters.get("z", []), dtype=float).flatten()
    has_ions = charges.size == int(mixture.ncomp) and bool(np.any(np.abs(charges) > 1.0e-12))
    if token in {"electrolyte_bubble_pressure", "electrolyte_bubble", "bubble_pressure"}:
        return {"route": "electrolyte_bubble", "reason": "requested electrolyte bubble-pressure path"}
    if token in {"electrolyte_lle", "electrolyte_lle_flash", "electrolyte_stability"}:
        return {"route": "electrolyte_lle", "reason": "requested electrolyte liquid-liquid path"}
    if token in {"lle", "lle_flash"}:
        return {"route": "neutral_lle", "reason": "requested neutral liquid-liquid path"}
    if token in {"tp_flash", "vle", "vle_flash"}:
        return {"route": "neutral_vle", "reason": "requested vapor-liquid path"}
    if token == "stability":
        return {"route": "neutral_tpd", "reason": "requested neutral stability path"}
    if token == "auto":
        if has_ions:
            return {"route": "electrolyte_lle", "reason": "ion-containing mixture"}
        if backend == "neutral_lle":
            return {"route": "neutral_lle", "reason": "neutral liquid-only backend requested"}
        return {"route": "neutral_vle", "reason": "neutral mixture default"}
    return {"route": "unsupported", "reason": "unsupported equilibrium kind"}
