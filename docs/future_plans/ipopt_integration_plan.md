# IPOPT Integration Plan for ePC-SAFT

## 1. Current ePC-SAFT State

The current package already has two distinct nonlinear workflows:

- Density solving is handled natively in the C++ core. The main entrypoint is [`den_cpp(...)`](C:/Users/Tanner/Documents/git/ePC-SAFT/src/epcsaft/epcsaft_electrolyte.cpp#L3899), which scans reduced-density space, refines root brackets, solves bracketed roots with Brent's method, and then validates candidate roots with pressure residual, positive \(\partial P / \partial \rho\), and residual Gibbs energy checks. The bound state path calls this through [`ePCSAFTStateNative::density()`](C:/Users/Tanner/Documents/git/ePC-SAFT/src/epcsaft/epcsaft_electrolyte.cpp#L4428).
- Regression is currently implemented in Python in [`src/epcsaft/regression.py`](C:/Users/Tanner/Documents/git/ePC-SAFT/src/epcsaft/regression.py#L717). The optimizer wrapper uses SciPy `least_squares(...)` for bounded nonlinear least-squares solves, and `fit_pure_neutral(...)` repeatedly evaluates the runtime through `mixture.state(...)` calls to build residuals from density and vapor-pressure consistency targets.
- Multiphase and LLE package workflows are intentionally not active in the current public package. Several analysis helpers now raise explicit `NotImplementedError` messages noting that the legacy multiphase workflow was removed and will need to be rewritten later in a new form.

The practical consequence is that ePC-SAFT does not currently need a general large-scale NLP solver for its primary runtime path. It does, however, have a natural future use case for one in regression and in later equilibrium formulations.

## 2. How IPOPT Is Built and Used

IPOPT is a large-scale nonlinear programming solver from COIN-OR. It is designed for smooth nonlinear programs with variable bounds, equality constraints, inequality constraints, sparse Jacobians, and optional sparse Hessians.

From the official IPOPT installation guidance:

- IPOPT requires BLAS and LAPACK.
- IPOPT also requires at least one sparse symmetric indefinite linear solver.
- On Windows, official build guidance is oriented around MSYS2/MinGW or COIN-OR CoinBrew-based builds.
- Optional linear-solver choices materially affect performance and robustness.

For code integration, the standard C++ pattern is through the `Ipopt::TNLP` interface plus an `IpoptApplication` instance. The problem class is expected to provide:

1. Problem dimensions.
2. Variable bounds.
3. Constraint bounds.
4. Initial values.
5. Jacobian sparsity structure and values.
6. Hessian sparsity structure and values, unless a quasi-Newton approximation is used.

In other words, IPOPT is not just a root finder. It expects a full nonlinear programming representation with explicit decision variables, constraint residuals, and first derivatives. This matters for ePC-SAFT because using IPOPT well generally means exposing an optimization model, not simply swapping it in for a scalar solver call.

IPOPT can also run without exact second derivatives. The official special-features documentation notes that setting `hessian_approximation=limited-memory` activates a limited-memory quasi-Newton approximation, which removes the need to provide an exact Hessian implementation. That is especially relevant for early ePC-SAFT integration because first derivatives are more realistic to add than full exact Hessians.

## 3. How Pyomo Integrates IPOPT

Pyomo currently supports IPOPT in two materially different ways.

### 3.1 File-based `ipopt` executable path

Pyomo has solver interfaces that write an `.nl` file and call an external `ipopt` executable. Relevant implementations include:

- [`pyomo.contrib.appsi.solvers.ipopt.Ipopt`](C:/Users/Tanner/Documents/git/pyomo/pyomo/contrib/appsi/solvers/ipopt.py#L137)
- [`pyomo.contrib.solver.solvers.ipopt.Ipopt`](C:/Users/Tanner/Documents/git/pyomo/pyomo/contrib/solver/solvers/ipopt.py#L237)

In this path, Pyomo:

- constructs a symbolic optimization model,
- writes the AMPL `.nl` representation,
- launches `ipopt`,
- parses the `.sol` output,
- and loads primal and dual results back into the model.

This route is appropriate when the entire model can be represented as standard Pyomo expressions or AMPL-compatible external functions.

### 3.2 In-process `cyipopt` path

Pyomo also supports in-process solves through `cyipopt` and PyNumero. Relevant code paths include:

- [`CyIpoptProblemInterface`](C:/Users/Tanner/Documents/git/pyomo/pyomo/contrib/pynumero/interfaces/cyipopt_interface.py#L66)
- [`CyIpoptNLP`](C:/Users/Tanner/Documents/git/pyomo/pyomo/contrib/pynumero/interfaces/cyipopt_interface.py#L255)
- [`PyomoCyIpoptSolver`](C:/Users/Tanner/Documents/git/pyomo/pyomo/contrib/pynumero/algorithms/solvers/cyipopt_solver.py#L245)

This route is more direct. The solver interface provides IPOPT with:

- variable values,
- bounds,
- constraints,
- gradients,
- Jacobians,
- and optionally Hessians

through Python callback methods instead of an external `.nl` file and subprocess.

### 3.3 External compiled-model path

For models that are not conveniently expressed as pure symbolic Pyomo equations, Pyomo provides the `ExternalGreyBoxModel` abstraction in [`external_grey_box.py`](C:/Users/Tanner/Documents/git/pyomo/pyomo/contrib/pynumero/interfaces/external_grey_box.py#L89). This is the most relevant Pyomo pattern for future ePC-SAFT equilibrium work.

The grey-box path is built for models that provide:

- named inputs,
- equality residuals and/or outputs,
- Jacobians of those residuals or outputs,
- and optionally Hessian products.

Pyomo’s own parameter-estimation example shows both the pure symbolic IPOPT route and the `cyipopt` external grey-box route side by side in [`perform_estimation.py`](C:/Users/Tanner/Documents/git/pyomo/pyomo/contrib/pynumero/examples/external_grey_box/param_est/perform_estimation.py#L58).

For ePC-SAFT, this matters because future flash, LLE, and chemical-equilibrium formulations will likely involve native EOS residuals that are easier to expose as grey-box callbacks than as a large pure-symbolic Pyomo expression tree.

## 4. Recommended ePC-SAFT Strategy

### 4.1 Do not replace the current density solver with IPOPT

The current density workflow should stay as-is.

Reasons:

- It is a one-dimensional root-finding problem, not a general NLP.
- The current implementation already contains tailored logic for scanning reduced density, bracketing roots, rejecting discontinuities, filtering metastable roots, and selecting phase-appropriate solutions.
- Replacing this with IPOPT would require reframing density as an optimization problem with artificial variables and constraints, which adds complexity without obvious benefit.
- The current path is native, local, and purpose-built for this specific thermodynamic closure problem.

So IPOPT should not be the default solver for `state(..., P=..., phase=...).density()`.

### 4.2 First use IPOPT as an optional regression backend

Regression is the right first target.

Reasons:

- The regression stack is already in Python.
- The current structure in [`_run_optimizer(...)`](C:/Users/Tanner/Documents/git/ePC-SAFT/src/epcsaft/regression.py#L717) cleanly isolates the optimizer call.
- The residual evaluation is already centralized and can be reused.
- IPOPT fits naturally as a bounded NLP backend for minimizing \(\frac{1}{2}\|r(x)\|_2^2\).

The first implementation should keep SciPy as the default and add IPOPT only as an optional backend.

### 4.3 Keep the first IPOPT implementation Python-level and optional

The first IPOPT implementation should use `cyipopt`, not a direct native C++ link against Ipopt.

Reasons:

- The current extension build is simple and lightweight: Cython plus the local C++ EOS core, with no IPOPT-specific link dependencies in [`pyproject.toml`](C:/Users/Tanner/Documents/git/ePC-SAFT/pyproject.toml#L22) or [`setup.py`](C:/Users/Tanner/Documents/git/ePC-SAFT/setup.py#L24).
- A direct Ipopt native link would make Windows packaging and developer setup materially harder.
- A Python-level `cyipopt` backend can remain optional and isolated from standard users.
- Using IPOPT with limited-memory Hessians is feasible even before exact Hessians are exposed.

### 4.4 Use Pyomo `ExternalGreyBoxModel` plus `cyipopt` later for equilibrium problems

For future flash, LLE, and chemical-equilibrium work, the best architectural fit is:

- Pyomo model for high-level decision variables and objective,
- `ExternalGreyBoxModel` for native EOS residual blocks,
- `cyipopt` for in-process nonlinear solves.

This is a much better fit than trying to route everything through the file-based `.nl` path, especially when EOS residuals and derivative callbacks will be computed in custom native code.

## 5. Implementation Roadmap

### Phase 1: optional `cyipopt` backend for `fit_pure_neutral(...)`

Add an optimizer selection layer so regression can choose between:

- SciPy `least_squares` as the default backend
- optional `cyipopt` as an advanced backend

Suggested implementation shape:

- keep the current residual-building logic unchanged,
- form the scalar objective \(\frac{1}{2}\|r(x)\|_2^2\),
- compute first derivatives with finite-difference residual Jacobians in the first version,
- configure IPOPT with `hessian_approximation=limited-memory`.

This phase should stay entirely in Python and should not modify the native C++ extension linkage.

### Phase 2: internal optimizer abstraction

After Phase 1 works, introduce a small internal abstraction layer so SciPy and IPOPT share the same:

- decision-variable normalization,
- bounds handling,
- multistart logic,
- residual assembly,
- result packaging.

This would reduce optimizer-specific branching in [`regression.py`](C:/Users/Tanner/Documents/git/ePC-SAFT/src/epcsaft/regression.py#L757) and make later nonlinear backends easier to add.

### Phase 3: experimental Pyomo grey-box equilibrium prototype

Once regression is stable, add an experimental equilibrium prototype using:

- Pyomo,
- `ExternalGreyBoxModel`,
- and `cyipopt`.

The first equilibrium target should be a square residual system such as a simple flash or binary LLE prototype with:

- phase compositions,
- density closures,
- mass-balance constraints,
- fugacity or electrochemical-potential equalities.

This should be explicitly experimental and kept separate from the current public runtime API.

### Phase 4: native C++ Ipopt linkage only if needed

Only consider direct native C++ Ipopt linkage if the Python-level route later proves too limited in performance, derivative fidelity, or solver-control requirements.

If this phase is ever reached, it should be treated as a separate packaging effort because it would require:

- new build detection logic,
- additional link libraries,
- Windows-specific toolchain handling,
- and likely optional native build flags.

## 6. Environment Notes

The current local `ePC-SAFT` environment was checked directly.

Observed status in the active `ePC-SAFT` Conda environment:

- SciPy is available.
- `pyomo` is not installed.
- `cyipopt` is not installed.
- Python module `ipopt` is not installed.
- no `ipopt` executable was found on `PATH`.

So IPOPT support should be documented as an optional setup path, not assumed to exist in the base developer environment.

For the first optional setup path, conda-forge is the most practical default because current `cyipopt` documentation presents Conda as the easiest cross-platform installation route. That should likely be the first documented path if and when optional IPOPT support is added to this repo.

## 7. References

- IPOPT installation and build documentation: <https://coin-or.github.io/Ipopt/INSTALL.html>
- IPOPT interface documentation: <https://coin-or.github.io/Ipopt/INTERFACES.html>
- IPOPT options documentation: <https://coin-or.github.io/Ipopt/OPTIONS.html>
- IPOPT special features, including limited-memory Hessians: <https://coin-or.github.io/Ipopt/SPECIALS.html>
- Pyomo `cyipopt` interface docs: <https://pyomo.readthedocs.io/en/6.9.3/api/pyomo.contrib.pynumero.interfaces.cyipopt_interface.html>
- Pyomo external grey-box model docs: <https://pyomo.readthedocs.io/en/latest/api/pyomo.contrib.pynumero.interfaces.external_grey_box.ExternalGreyBoxModel.html>
- `cyipopt` installation docs: <https://cyipopt.readthedocs.io/en/stable/install.html>
