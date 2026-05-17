from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

TEXT_SUFFIXES = {
    ".bib",
    ".cfg",
    ".cmake",
    ".cpp",
    ".h",
    ".hpp",
    ".ini",
    ".md",
    ".py",
    ".pyi",
    ".rst",
    ".sh",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
    ".ps1",
}

EXCLUDED_PREFIXES = (
    "docs/papers/",
)
EXCLUDED_PATHS: set[str] = set()


def _tracked_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    paths: list[Path] = []
    for line in result.stdout.splitlines():
        rel = line.strip().replace("\\", "/")
        if not rel:
            continue
        if rel in EXCLUDED_PATHS or any(rel.startswith(prefix) for prefix in EXCLUDED_PREFIXES):
            continue
        path = REPO_ROOT / rel
        if not path.exists():
            continue
        if path.suffix.lower() in TEXT_SUFFIXES:
            paths.append(path)
    return paths


def _blocked_terms() -> tuple[str, ...]:
    first = "backend" + "_" + "unavailable"
    second = "finite" + "_" + "difference"
    third = "finite" + "-" + "difference"
    fourth = "finite" + " " + "difference"
    fifth = "central" + "_" + "perturbation"
    sixth = "source" + "_" + "perturbation"
    seventh = "perturbation" + "_" + "derivative"
    eighth = "numerical" + "diff"
    ninth = "levenberg" + "_" + "marquardt"
    tenth = "manual" + " numeric " + "perturbation"
    eleventh = "numerical" + " perturbation"
    twelfth = "numeric" + " perturbation"
    thirteenth = "perturbation" + " jacobian"
    fourteenth = "perturbation" + " derivative"
    fifteenth = "perturbation" + "-" + "jacobian"
    sixteenth = "perturbation" + "-" + "derivative"
    seventeenth = "perturbation" + "-" + "derived"
    eighteenth = "numerical" + "_" + "derivative"
    nineteenth = "numerical" + "_" + "jacobian"
    twentieth = "numerical" + " derivative"
    twenty_first = "numeric" + " derivative"
    twenty_second = "numerical" + " jacobian"
    twenty_third = "numeric" + " jacobian"
    twenty_fourth = "perturbation" + " differencing"
    twenty_fifth = "differential_mode" + '": "' + "numerical"
    twenty_sixth = "residual" + "_evaluation" + "_only"
    twenty_seventh = "central" + " perturbation"
    twenty_eighth = "central" + " perturbations"
    twenty_ninth = "density" + "-perturbation"
    thirtieth = "local density " + "perturbation"
    thirty_first = "perturbation" + " slope"
    thirty_second = "perturbation" + " slopes"
    thirty_third = "numeric" + "_diff"
    thirty_fourth = "numeric" + "-" + "diff"
    thirty_fifth = "numeric" + "diff"
    thirty_sixth = "numerical" + "_diff"
    thirty_seventh = "numerical" + "-" + "diff"
    thirty_eighth = "numeric" + " differencing"
    thirty_ninth = "numerical" + " differencing"
    fortieth = "numeric" + " differentiation"
    forty_first = "numerical" + " differentiation"
    forty_second = "numerical" + "-" + "differentiation"
    forty_third = "initial" + "_phases" + "_from" + "_result"
    forty_fourth = "charge" + "_neutral" + "_lle" + "_seed"
    forty_fifth = "solvent" + "_endpoint" + "_seed"
    forty_sixth = "electrolyte" + "_seeds"
    forty_seventh = "route" + "_pending"
    forty_eighth = "route" + "-" + "pending"
    forty_ninth = "route" + " pending"
    fiftieth = "pending" + " route"
    fifty_first = "pending" + "-" + "route"
    return (
        first,
        second,
        third,
        fourth,
        fifth,
        sixth,
        seventh,
        eighth,
        ninth,
        tenth,
        eleventh,
        twelfth,
        thirteenth,
        fourteenth,
        fifteenth,
        sixteenth,
        seventeenth,
        eighteenth,
        nineteenth,
        twentieth,
        twenty_first,
        twenty_second,
        twenty_third,
        twenty_fourth,
        twenty_fifth,
        twenty_sixth,
        twenty_seventh,
        twenty_eighth,
        twenty_ninth,
        thirtieth,
        thirty_first,
        thirty_second,
        thirty_third,
        thirty_fourth,
        thirty_fifth,
        thirty_sixth,
        thirty_seventh,
        thirty_eighth,
        thirty_ninth,
        fortieth,
        forty_first,
        forty_second,
        forty_third,
        forty_fourth,
        forty_fifth,
        forty_sixth,
        forty_seventh,
        forty_eighth,
        forty_ninth,
        fiftieth,
        fifty_first,
    )


def _source_blocked_terms(rel: str) -> tuple[str, ...]:
    if not rel.startswith(("src/", "tests/", "scripts/")):
        return ()
    terms = (
        "fall" + "back",
        "return" + "_best" + "_effort",
        "best" + "_effort",
        "_solve" + "_electrolyte" + "_bubble" + "_native",
        "native" + "_log" + "_pressure" + "_bisection",
        "hessian" + "_backend",
        "multi" + "start",
        "line" + "_search",
        "soft" + "_start",
        "damp" + "ing",
        "native" + "_scalar" + "_binary" + "_activity",
        "binary" + "_log" + "_amounts",
        "derivative" + "_status",
        "least" + "_squares",
        "differential" + "_evolution",
        "minimize" + "_scalar",
        "new" + "ton",
        "root" + "_scalar",
        "brent" + "q",
        "brent" + "h",
        "without" + " optimizer",
        "golden" + "_section",
        "golden" + "-" + "section",
        "evaluate" + "_generic" + "_regression" + "_derivatives",
        "residual" + "_score" + "_native",
        "sci" + "py.optimize",
        "np.linalg." + "lstsq",
        "numpy.linalg." + "lstsq",
        "implicit" + "_sensitivity" + "_status",
        "not" + " implemented",
        "association" + "_solver" + "_status",
        "ceres support" + " is not enabled",
        "ceres and cppad support" + " are not enabled",
        "ceres" + "_disabled",
        "cppad" + "_disabled",
        "cppad support" + " is disabled",
        "optional" + "_dependencies",
        "dependency" + "_not" + "_compiled",
        "reactive" + "_flash" + "_tp",
        "_normalize" + "_neutral" + "_initial" + "_phases",
        "_electrolyte" + "_initial" + "_phase" + "_seed",
        "material" + "_balanced" + "_initial" + "_phases",
        "solver" + "_seed" + "_name",
        "optimizer" + "_backend" + "\": \"not" + "_applicable",
        "optimizer" + "_backend == \"not" + "_applicable",
        "derivative" + "_backend" + "\": \"not" + "_applicable",
        "derivative" + "_backend == \"not" + "_applicable",
        "jacobian" + "_backend" + "\": \"not" + "_applicable",
        "jacobian" + "_backend == \"not" + "_applicable",
        "jacobian" + "_backend = \"not" + "_applicable",
        "ascani" + "_benchmark" + "_attempt",
        "not" + "_applicable" + "_to" + "_neutral" + "_route",
        "requires" + "_ipopt" + "_build",
        "unsupported" + "_reaction" + "_scopes",
        "not" + "_configured",
        "native" + "_ipopt" + "_equilibrium" + "_nlp" + "_required",
        "not" + "_available",
        "validated" + "_for" + "_pending" + "_ipopt" + "_route",
        "route" + "_gated" + "_standard" + "_states",
        "explicit" + "_cppad" + "_request" + "_raises" + "_until" + "_implemented",
        "allow" + "_unsupported" + "_parameters",
        "allow" + "_without" + "_direct" + "_data",
        "partially" + "_available",
        "not" + "_applicable",
        "initial" + "_y" + "_vap",
        "options.initial" + "_pressure",
        "electrolytebubbleoptions(initial" + "_pressure",
        "_with" + "_bubble" + "_continuation",
        "initial" + "_phases" + "=none",
        "initial" + "_phases" + " = none",
        "phase" + "_kwargs.get(\"initial" + "_phases\")",
        "phase" + "_kwargs={\"initial" + "_phases\"",
        "package" + "-owned alternate",
        "from" + "_legacy" + "_record",
        "p" + "_seed",
        "density" + "_seed" + "_aliases",
        "stability" + "_precheck",
        "include" + "_phase" + "_diagnostics",
        "python" + "_sweep",
        "selected" + "_solver" + "_backend" + "\": \"not" + "_run",
        "residual" + "_surface" + "_only",
        "hessian" + "_strategy",
        "hessian" + "_strategies",
        "hessian" + "_kind",
        "exact" + "_hessian" + "_required",
        "exact" + "_hessian" + "_available",
        "hessian" + "_callback" + "_available",
        "hessian" + "_includes" + "_second" + "_residual" + "_derivatives",
        "requires" + "_exact" + "_hessian",
        "solver" + "_hessian" + "_mode",
        "activity" + "_derivative" + "_policy",
        "activity" + "_fixed" + "_point",
        "activity" + "_or" + "_fugacity" + "_terms" + "_in" + "_residual",
        "activity" + "_derivative" + "_in" + "_jacobian",
        "continuation" + "_state",
        "continuation" + "_used",
        "covariance" + "_available" + "\": false",
        "fit" + "_success" + "\": none,",
        "fit" + "_message" + "\": \"\",",
        "termination" + "_reason" + "\": \"objective" + "_only\",",
        "fit" + "_iterations" + "\": none,",
        "objective" + "_initial" + "\": none,",
        "gradient" + "_norm" + "\": none,",
        "step" + "_norm" + "\": none,",
        "diagnostic" + "_residual" + "_score",
        "diagnostic" + "_fields",
        "package" + "_status",
        "comparison" + "_status",
        "package" + "_missing",
        "diagnostic" + "_failure",
        "best" + "_failure" + "_reason",
        "native" + "_status",
        "derivative" + "_gap" + "_status",
        "\"warm" + "_start" + "_source\"",
        "objective" + "_seed" + "_hits",
        "row" + "_seed" + "_hits",
        "baseline" + "_repeat",
        "baseline" + "_warmup",
        "production" + "_eos" + "_coverage",
        "unsupported" + "_derivative" + "_behavior",
        "full" + "_constrained" + "_nlp" + "_available",
        "fit" + "_reactive" + "_electrolyte" + "_parameters",
        "reactive" + "regression" + "fit" + "result",
        "unsupported" + " reaction constant convention",
        "reaction standard state contains an " + "unsupported code",
    )
    if rel.startswith(("src/", "scripts/")):
        terms += (
            "warm" + "_start" + "_rows",
            "warm" + "_start" + "_objective",
            "warm" + "_start" + "_used",
            "warm" + "_start" + "_failed",
            "objective" + "_seed",
            "row" + "_seed",
        )
    if rel.startswith("src/"):
        terms += ("fit" + "_route",)
        terms += ("activity" + "_coefficient" + "_term" + "_decomposition" + "_available",)
        terms += ("not" + "_validated" + "_for" + "_production",)
    if rel in {"src/epcsaft/reactive_regression.py", "scripts/benchmarks/helpers/reactive_regression.py"}:
        terms += (
            "warm" + "_start" + "_source",
            "context" + "_cache",
        )
    return terms


def main() -> int:
    terms = _blocked_terms()
    matches: list[str] = []
    for path in _tracked_files():
        rel = path.relative_to(REPO_ROOT).as_posix()
        text = path.read_text(encoding="utf-8", errors="ignore").lower()
        path_terms = terms + _source_blocked_terms(rel)
        for line_number, line in enumerate(text.splitlines(), start=1):
            if any(term in line for term in path_terms):
                matches.append(f"{rel}:{line_number}: blocked solver/derivative gate term")

    if matches:
        print("Strict solver derivative text gate failed:", file=sys.stderr)
        print("\n".join(matches), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
