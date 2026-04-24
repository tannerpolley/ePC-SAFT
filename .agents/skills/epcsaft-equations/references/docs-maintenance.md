# Docs And Skill Maintenance

## What Belongs Here

Keep detailed ePC-SAFT science, equation ownership, and equation-to-code navigation in this skill. Keep AGENTS.md focused on operating rules, workflow commands, delegation policy, and memory/commit policy.

Update this skill when any of these change:

- The public Python to native path changes, especially `src/epcsaft/epcsaft.py`, `src/epcsaft/bindings.cpp`, or `epcsaft._core`.
- Native equation ownership moves between `src/epcsaft/native/**` files.
- Contribution families, payload keys, EqID ownership comments, or density/phase closure behavior change.
- `docs/latex/equations.tex`, `docs/equations.md`, or `docs/equations_registry.yaml` change in a way that affects how agents should trace equations.
- Build/tooling changes affect equation-facing work, such as pybind/CMake source ownership or generated equation sync commands.

## What Belongs In User-Facing Docs

Update README/docs when users or developers need different instructions:

- Install/build/test commands change.
- Public API behavior changes.
- Parameter dataset, regression, or user-options workflows change.
- A documented equation or derived property behavior changes.

Prefer concise current-state docs. Remove obsolete workaround history unless it is still needed to prevent a recurring mistake.

## Minimal Maintenance Checklist

For equation or native-runtime changes:

1. Check whether `docs/latex/equations.tex` is still the source of truth for the changed behavior.
2. Check whether `docs/equations.md` and `docs/equations_registry.yaml` need regeneration.
3. Check whether `references/code-map.md` still points to the correct owner files.
4. Check whether README or Sphinx docs need user-facing updates.
5. Record only durable workflow or architecture changes in `docs/.codex-journal/project_memory.md`.

## Native Build Coordination

Equation and native-runtime validation may import `epcsaft._core`, which can lock the generated `_core*.pyd` on Windows. Prefer normal native builds for routine checks, and avoid clean/repair rebuilds while tests, Python REPLs, IDE run configurations, or sub-agents may have the extension loaded. If a clean repair is needed, coordinate it through the main orchestrator first.
