from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import _common as common
from scripts._env import require_epcsaft_install

require_epcsaft_install()

from epcsaft.parameters import get_prop_dict

DATA_PATH = Path(__file__).with_name("1-butanol-NH4Cl-water-LLE.csv")
SPECIES = ["H2O", "Butanol", "NH4+", "Cl-"]
MW = np.asarray([18.0153e-3, 74.1216e-3, 18.038e-3, 35.453e-3], dtype=float)
MW_NH4CL = float(MW[2] + MW[3])
IDX = {name: i for i, name in enumerate(SPECIES)}


def _split_nh4cl_mass(total_salt_mass_fraction: float) -> tuple[float, float]:
    total = float(total_salt_mass_fraction)
    return total * float(MW[2] / MW_NH4CL), total * float(MW[3] / MW_NH4CL)


def _mass_to_mole_fraction(w: np.ndarray) -> np.ndarray:
    n = np.asarray(w, dtype=float) / MW
    return n / np.sum(n)


def _mole_to_mass_fraction(x: np.ndarray) -> np.ndarray:
    w = np.asarray(x, dtype=float) * MW
    return w / np.sum(w)


def _solve_lle(feed_mass_fraction: np.ndarray) -> dict | None:
    raise NotImplementedError("The legacy multiphase LLE workflow has been removed and will be rewritten later.")


@lru_cache(maxsize=1)
def load_experimental_rows() -> tuple[dict, ...]:
    rows: list[dict] = []
    with DATA_PATH.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                w_nh4cl_aq = float(row["w%_NH4Cl in aqueous phase"]) / 100.0
                w_buoh_org = float(row["w%_1-butanol in 1-butanol phase"]) / 100.0
                w_buoh_aq = float(row["w%_1-butanol in aqueous phase"]) / 100.0
            except (KeyError, TypeError, ValueError):
                continue

            m_nh4, m_cl = _split_nh4cl_mass(w_nh4cl_aq)
            aqueous = np.asarray(
                [
                    1.0 - w_nh4cl_aq - w_buoh_aq,
                    w_buoh_aq,
                    m_nh4,
                    m_cl,
                ],
                dtype=float,
            )
            # The provided Figure 6 data only reports NH4Cl in the aqueous phase and 1-butanol in both phases.
            # For plotting the experimental ternary points, the remaining organic fraction is treated as water.
            organic = np.asarray([1.0 - w_buoh_org, w_buoh_org, 0.0, 0.0], dtype=float)

            rows.append(
                {
                    "w_nh4cl_aq": w_nh4cl_aq,
                    "w_buoh_org": w_buoh_org,
                    "w_buoh_aq": w_buoh_aq,
                    "aqueous_mass_fraction": aqueous,
                    "organic_mass_fraction": organic,
                }
            )
    return tuple(rows)


def _objective(target: dict, solved: dict) -> float:
    aq = np.asarray(solved["aqueous_mass_fraction"], dtype=float)
    org = np.asarray(solved["organic_mass_fraction"], dtype=float)
    target_aq = np.asarray(target["aqueous_mass_fraction"], dtype=float)
    target_org = np.asarray(target["organic_mass_fraction"], dtype=float)

    aq_salt = float(aq[IDX["NH4+"]] + aq[IDX["Cl-"]])
    org_salt = float(org[IDX["NH4+"]] + org[IDX["Cl-"]])
    err = 0.0
    err += 3.0 * (aq_salt - float(target["w_nh4cl_aq"])) ** 2
    err += 1.5 * (float(aq[IDX["Butanol"]]) - float(target_aq[IDX["Butanol"]])) ** 2
    err += 0.8 * (float(org[IDX["Butanol"]]) - float(target_org[IDX["Butanol"]])) ** 2
    err += 0.2 * org_salt**2
    if float(org[IDX["Butanol"]]) < float(aq[IDX["Butanol"]]):
        err += 10.0
    return float(err)


@lru_cache(maxsize=1)
def solve_model_rows() -> tuple[dict, ...]:
    solved_rows: list[dict] = []
    for idx, row in enumerate(load_experimental_rows()):
        aq = np.asarray(row["aqueous_mass_fraction"], dtype=float)
        org = np.asarray(row["organic_mass_fraction"], dtype=float)

        # The LLE solver in this package requires ions, so use a trace NH4Cl amount for the salt-free anchor.
        if float(row["w_nh4cl_aq"]) <= 1.0e-12:
            trace_salt = 1.0e-6
            m_nh4, m_cl = _split_nh4cl_mass(trace_salt)
            aq_eval = aq.copy()
            aq_eval[IDX["H2O"]] -= trace_salt
            aq_eval[IDX["NH4+"]] += m_nh4
            aq_eval[IDX["Cl-"]] += m_cl
        else:
            aq_eval = aq

        best: tuple[float, float, dict] | None = None
        for lam in (0.85, 0.90, 0.95, 0.98, 0.995):
            feed = lam * aq_eval + (1.0 - lam) * org
            solved = _solve_lle(feed)
            if solved is None or solved["aqueous_mass_fraction"] is None or solved["organic_mass_fraction"] is None:
                continue
            score = _objective(row, solved)
            if best is None or score < best[0]:
                best = (score, float(lam), solved)

        if best is None:
            continue

        score, lam, solved = best
        aq_model = np.asarray(solved["aqueous_mass_fraction"], dtype=float)
        org_model = np.asarray(solved["organic_mass_fraction"], dtype=float)
        solved_rows.append(
            {
                "row_id": idx,
                "feed_lambda": lam,
                "objective": score,
                "aqueous_mass_fraction": aq_model,
                "organic_mass_fraction": org_model,
                "w_nh4cl_aq": float(aq_model[IDX["NH4+"]] + aq_model[IDX["Cl-"]]),
                "w_buoh_org": float(org_model[IDX["Butanol"]]),
                "w_buoh_aq": float(aq_model[IDX["Butanol"]]),
                "beta_organic": float(solved["beta_organic"]),
                "beta_aqueous": float(solved["beta_aqueous"]),
            }
        )

    solved_rows.sort(key=lambda item: item["w_nh4cl_aq"])
    return tuple(solved_rows)
