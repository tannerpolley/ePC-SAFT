# ePC-SAFT Codex Environment

This ignored local folder intentionally has one environment: `environment.toml`.

The setup path is:

```powershell
uv sync --no-install-project
pwsh.exe -NoProfile -ExecutionPolicy Bypass -File .codex/environments/setup.ps1
```

The native build action intentionally uses the script default ``--profile fast``:
Ceres, CppAD, and native Ipopt are enabled when available. The setup helper resolves
``EPCSAFT_IPOPT_ROOT`` first, then ``EPCSAFT_PEP517_IPOPT_ROOT``, then the local
Windows SDK default ``%USERPROFILE%\Documents\deps\ipopt-msvc``. It also adds
the Ipopt ``bin`` directory to ``PATH`` and ``EPCSAFT_RUNTIME_DLL_DIRS`` so
``doctor.py --require-ipopt`` can load the native extension in a fresh Codex
worktree. Use ``uv run python scripts/dev/validate_project.py ceres-cppad`` when
the task needs the focused Ceres regression/backend slice.
Native Ipopt discovery outside this Windows setup remains explicit
system-dependency work. Use ``EPCSAFT_IPOPT_ROOT=<Ipopt-root>`` or
``EPCSAFT_PEP517_IPOPT_ROOT=<Ipopt-root>`` when the local default is not
present. Use ``build_epcsaft.py --ipopt-dir <IpoptConfig-dir>`` manually when
only an ``IpoptConfig.cmake`` directory is available.
Existing ``build/dev`` trees configured with MinGW/Strawberry cannot be switched
in place to the MSVC Ipopt toolchain; use a fresh Codex worktree or run the
coordinated clean build repair before enabling this setup on an old checkout.
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
