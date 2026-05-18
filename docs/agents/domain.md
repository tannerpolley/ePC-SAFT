# Domain Docs

How engineering skills should consume this repo's domain documentation when exploring the codebase.

## Before exploring, read these

- **`CONTEXT.md`** at the repo root.
- **`docs/adr/`** at the repo root. Start with `docs/adr/README.md`, then read any numbered ADRs relevant to the current topic.
- **`CONTEXT-MAP.md`** only if this repo later becomes multi-context and adds one.

If a needed concept is missing from `CONTEXT.md`, note the gap or use `grill-with-docs` to resolve the terminology before naming a new package concept in durable output.

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
