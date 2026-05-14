# ePC-SAFT Codex Environment

This ignored local folder intentionally has one environment: `environment.toml`.

The setup path is:

```powershell
uv sync --no-install-project
uv run python scripts/dev/build_epcsaft.py
uv run python scripts/dev/doctor.py
```

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
