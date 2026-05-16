# ePC-SAFT Codex Environment

This ignored local folder intentionally has one environment: `environment.toml`.

The setup path is:

```powershell
uv sync --no-install-project
uv run python scripts/dev/build_epcsaft.py
uv run python scripts/dev/doctor.py
```

The native build action intentionally uses the script default ``--profile fast``:
Ceres and CppAD are both required, while Ipopt is explicit adapter work.
Use ``uv run python scripts/dev/validate_project.py ceres-cppad`` when the task
needs the focused Ceres regression/backend slice.
Native Ipopt discovery is explicit system-dependency work. Use
``uv run python scripts/dev/build_epcsaft.py --profile ipopt --ipopt-root <Ipopt-root>``
or ``uv run python scripts/dev/build_epcsaft.py --profile ipopt --ipopt-dir <IpoptConfig-dir>``
only for native Ipopt adapter development or validation.
This local dev-script default matches package installs for required Ceres/CppAD
dependencies.
For repeated full package installs, build reusable Ceres once with
``uv run python scripts/dev/build_system_ceres.py --parallel 4`` and set the
printed ``EPCSAFT_PEP517_CERES_DIR`` / ``EPCSAFT_PEP517_BUILD_DIR`` variables.

The action list should stay lean and limited to the normal project workflow:

- `Sync Environment`
- `Doctor`
- `Build Native Extension`
- `Validate Quick`
- `Validate Confidence`
- `Build Docs`
- `Build Distribution`

Advanced commands such as fast native rebuilds, clean repair, profiling, equation registry maintenance, and targeted tests stay in `docs/pages/development_workflows.rst` and should be run manually when needed.

Maintenance rules:

- Keep this file and `environment.toml` current whenever script names, test paths, docs commands, or validation expectations change.
- Remove stale actions immediately when a workflow is deleted.
- Do not re-add plot gallery, gallery server, manifest, index, `--plots`, or obsolete test-slice actions to this package.
- Before handoff after workflow edits, parse `environment.toml` and search for removed scripts, removed flags, old test paths, and obsolete docs/gallery references.
