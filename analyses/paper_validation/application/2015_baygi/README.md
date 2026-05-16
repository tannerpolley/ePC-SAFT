# 2015 Baygi Analysis

Pure MEA validation and regression smoke workflow.

The figure workflow renders from tracked cached diagnostics by default. Live
pure-MEA saturation recomputation is route-gated until the package has a native
Ipopt bubble/dew implementation; no analysis-local saturation solve loop is
retained.

Entry point:

```powershell
uv run python analyses\paper_validation\application\2015_baygi\scripts\run_all.py
```
