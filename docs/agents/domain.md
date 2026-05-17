# Domain Docs

How engineering skills should consume this repo's domain documentation when exploring the codebase.

## Before exploring, read these

- **`CONTEXT.md`** at the repo root, when present.
- **`docs/adr/`** at the repo root, when present.
- **`CONTEXT-MAP.md`** only if this repo later becomes multi-context and adds one.

If any of these files do not exist, proceed silently. Do not flag their absence or suggest creating them upfront. The producer skill (`grill-with-docs`) creates them when terminology or architectural decisions are actually resolved.

## Layout

This is treated as a single-context repo:

```text
/
|-- CONTEXT.md
|-- docs/
|   `-- adr/
`-- src/
```

## Use the glossary's vocabulary

When output names a domain concept in an issue title, refactor proposal, hypothesis, or test name, use the term as defined in `CONTEXT.md`.

If the needed concept is not in the glossary yet, either reconsider whether the language belongs in this repo or note the gap for `grill-with-docs`.

## Flag ADR conflicts

If output contradicts an existing ADR, surface it explicitly rather than silently overriding it.
