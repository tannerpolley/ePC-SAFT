# ePC-SAFT Codex Environment

This ignored local folder intentionally has one environment: `environment.toml`.

The setup path is:

```powershell
uv sync --no-install-project
uv run python scripts/dev/build_epcsaft.py
uv run python scripts/dev/doctor.py
```

The native build action intentionally uses the script default ``--profile fast``
(Ceres OFF, CppAD ON). Ceres builds are explicit validation work, not setup.
Use ``uv run python scripts/dev/build_epcsaft.py --profile full`` or
``uv run python scripts/dev/validate_project.py ceres-cppad`` only when the task
needs Ceres regression/backend coverage.
This local dev-script default does not change package installs; editable, wheel,
and downstream path installs still inherit the CMake default with Ceres ON.

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
