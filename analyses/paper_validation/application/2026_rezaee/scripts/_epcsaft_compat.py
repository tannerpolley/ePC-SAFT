from __future__ import annotations

from typing import Any

import numpy as np


def _as_array(key: str, value: Any) -> Any:
    if key == "assoc_scheme":
        return value
    if isinstance(value, list):
        return np.asarray(value, dtype=float)
    return value


def _modern_elec_model(model: Any) -> Any:
    if not isinstance(model, dict):
        return model
    if (
        "rel_perm" in model
        or "include_born_model" in model
        or isinstance(model.get("DH_model"), dict)
        or isinstance(model.get("born_model"), dict)
    ):
        allowed = {"rel_perm", "hc_model", "disp_model", "assoc_model", "DH_model", "include_born_model", "born_model"}
        return {key: value for key, value in model.items() if key in allowed}

    born_options = model.get("born_diff_options", {})
    if not isinstance(born_options, dict):
        born_options = {}
    differential_mode = model.get("dielc_diff_mode", "analytical")
    born_differential_mode = model.get("born_diff_model", "analytical")
    return {
        "rel_perm": {
            "rule": model.get("dielc_rule", 1),
            "differential_mode": differential_mode,
        },
        "hc_model": {"dadx_differential_mode": "analytical"},
        "disp_model": {"dadx_differential_mode": "analytical"},
        "assoc_model": {"dadx_differential_mode": "analytical"},
        "DH_model": {
            "d_ion_mode": 1,
            "bjeruum_treatment": bool(model.get("bjeruum_treatment", False)),
            "mu_DH_model": {
                "differential_mode": differential_mode,
                "comp_dep_rel_perm": True,
                "include_sum_term": True,
            },
        },
        "include_born_model": bool(model.get("born_contrib", True)),
        "born_model": {
            "d_Born_mode": model.get("d_Born_mode", 0),
            "solvation_shell_model": bool(model.get("ssm_ds", False)),
            "dielectric_saturation": bool(model.get("dielectric_saturation", False)),
            "bulk_mode": model.get("eps_r_bulk", "mix"),
            "mu_born_model": {
                "differential_mode": born_differential_mode,
                "comp_dep_rel_perm": bool(born_options.get("include_dielc_conc_dep", True)),
                "include_sum_term": bool(born_options.get("include_sum_term", True)),
                "comp_dep_delta_d": bool(born_options.get("include_delta_d_i_conc_dep", False)),
            },
        },
    }


def _normalized_params(params: dict[str, Any]) -> dict[str, Any]:
    removed_keys = {
        "dipm",
        "dip_num",
        "bjeruum_treatment",
        "born_model",
        "born_diff_mode",
        "born_eps_mode",
        "DH_model",
        "dielc_rule",
        "dielc_diff_mode",
        "debug",
        "elec_model_preset",
    }
    out = {key: _as_array(key, value) for key, value in dict(params).items() if key not in removed_keys}
    if "elec_model" in out:
        out["elec_model"] = _modern_elec_model(out["elec_model"])
    return out
