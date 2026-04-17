# IPOPT Integration Plan for ePC-SAFT

## Status

This note is superseded by [docs/native-ipopt-regression-plan.md](C:/Users/Tanner/Documents/git/ePC-SAFT/docs/native-ipopt-regression-plan.md).

The regression cutover no longer follows the older optional Python `cyipopt` direction. The implemented package path is now:

- native C++ IPOPT ownership inside the main `epcsaft.epcsaft` extension
- public Python surface kept at `fit_pure_neutral(...)`, `FitProblem`, `FitResult`, `load_regression_records(...)`, and `write_fit_result(...)`
- IPOPT treated as a required build dependency for the package

## Current v1 Scope

The shipped first phase remains intentionally narrow:

- one neutral component only
- fitted parameters limited to \(m\), \(s\), and \(e\)
- weighted least-squares objective over liquid-density and pure-VLE fugacity-balance records
- exact first derivatives
- density retained as a nested native solve

The canonical architecture, build, and validation notes now live in [docs/native-ipopt-regression-plan.md](C:/Users/Tanner/Documents/git/ePC-SAFT/docs/native-ipopt-regression-plan.md).
