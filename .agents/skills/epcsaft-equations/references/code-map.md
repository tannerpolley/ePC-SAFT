# Code Map

## Python And Pybind Surface

- `src/epcsaft/__init__.py`: public package exports.
- `src/epcsaft/parameters.py`: dataset/user-option normalization, association parsing, ion/Born/dielectric options, and parameter payload construction.
- `src/epcsaft/epcsaft.py`: public native-backed classes, Python method behavior, contribution payload formatting, validation, and calls into `epcsaft._core`.
- `src/epcsaft/regression.py`: public regression helpers and native least-squares result normalization.
- `src/epcsaft/bindings.cpp`: pybind11 binding layer that exposes the private extension module `epcsaft._core`.

Public entry points to trace first:

- `ePCSAFTMixture.from_params(...)`
- `ePCSAFTMixture.from_dataset(...)`
- `ePCSAFTMixture.state(...)`
- `ePCSAFTState.residual_helmholtz(...)`
- `ePCSAFTState.compressibility_factor(...)`
- `ePCSAFTState.fugacity_coefficient(...)`
- `ePCSAFTState.activity_coefficient(...)`

## Native C++ Owners

- `src/epcsaft/native/epcsaft_parameter_setup.cpp`: shared parameter and mixture setup helpers, including segment/ion/Born radii and pair setup rules.
- `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp`: hard-chain state and derivative contribution logic.
- `src/epcsaft/native/contributions/epcsaft_contrib_disp.cpp`: dispersion contribution and derivative logic.
- `src/epcsaft/native/contributions/epcsaft_contrib_assoc.cpp`: association intermediate state, site fractions, and association contribution logic.
- `src/epcsaft/native/contributions/epcsaft_contrib_ion.cpp`: electrolyte / Debye-Huckel contribution logic.
- `src/epcsaft/native/contributions/epcsaft_contrib_born.cpp`: Born contribution and solvation-shell / dielectric-saturation support.
- `src/epcsaft/native/contributions/epcsaft_contrib_internal.h`: shared contribution result/intermediate structs and contribution function declarations.
- `src/epcsaft/native/epcsaft_ares.cpp`: residual Helmholtz aggregation and derivative ownership.
- `src/epcsaft/native/epcsaft_Z.cpp`: compressibility factor and pressure-related ownership.
- `src/epcsaft/native/epcsaft_density.cpp`: density closure and pressure-based state solve behavior.
- `src/epcsaft/native/epcsaft_mu.cpp`: residual chemical potential ownership.
- `src/epcsaft/native/epcsaft_fugcoef.cpp`: fugacity coefficient ownership.
- `src/epcsaft/native/epcsaft_activity.cpp`: component and mean-ionic activity coefficient ownership.
- `src/epcsaft/native/epcsaft_state.cpp`: native mixture/state lifecycle, cache behavior, and state property methods.
- `src/epcsaft/native/epcsaft_thermo.cpp`: residual thermodynamic property helpers.

## Mapping Pattern

For a property bug or equation change:

1. Identify the public Python method in `src/epcsaft/epcsaft.py` or `src/epcsaft/regression.py`.
2. Find the matching pybind binding in `src/epcsaft/bindings.cpp`.
3. Trace to the native owner file above.
4. If the behavior is contribution-specific, trace into `contributions/*` through `epcsaft_contrib_internal.h`.
5. Check whether the equation has an `EqID` in `docs/latex/equations.tex` and whether the generated registry points to the same owner.
