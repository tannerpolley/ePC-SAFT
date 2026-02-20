from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> None:
    base = Path(__file__).resolve().parent
    csv_path = base / "water-methanol.csv"
    out_path = base / "gibbs_transfer_K_Cl_water_methanol.png"

    x_vals: list[float] = []
    k_vals: list[float] = []
    cl_vals: list[float] = []

    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            x_vals.append(float(row["x_MeOH"]))
            k_vals.append(float(row["K+"]))
            cl_vals.append(float(row["Cl-"]))

    plt.figure(figsize=(8, 5))
    plt.plot(x_vals, k_vals, marker="o", linewidth=2, label="K+")
    plt.plot(x_vals, cl_vals, marker="s", linewidth=2, label="Cl-")
    plt.xlabel("Methanol mole fraction, x_MeOH")
    plt.ylabel("Gibbs transfer")
    plt.title("Gibbs Transfer in Water-Methanol: K+ and Cl-")
    plt.xlim(0.0, 1.0)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
