# T999: Audit After Manifest Slice

Task: `T999`
Kind: `judge`
Status: `current`

## Summary

The manifest/CLI/test slice and the docs discoverability slice are both valid and verified, but the full issue outcome is not complete yet. The package now has a truthful literature-suite inventory surface plus a maintainer-facing docs entrypoint, but the suite still has multiple follow-up benchmark anchors classified as blockers and one visible package-scope violation in the MEA-specific public helper.

## Details

- Verified results:
  - `tests/workflows/test_benchmark_literature_suite.py` passed.
  - The issue-listed regression checks still passed.
  - `uv run python scripts/validate_project.py docs` passed.
  - `uv run python scripts/benchmark_literature_suite.py --case figiel_2025_ssm_ds_born` rendered the new suite surface successfully.
- Why the goal is not complete:
  - several issue-scope anchors remain classified as blockers rather than implemented;
  - the MEA benchmark path still depends on an application-specific compatibility helper.
- Next safe task:
  - map the smallest generic-regression migration path that removes `fit_mea_co2_h2o_electrolyte` from the literature benchmark path without introducing a new application-specific API.

## Board Receipt Snippet

```yaml
receipt:
  result: not_complete
  note: notes/T999-audit-after-manifest-slice.md
  summary: "The manifest and docs slices are verified, but the goal is not complete: several literature anchors remain blocked and the MEA benchmark path still depends on an application-specific public helper."
```
