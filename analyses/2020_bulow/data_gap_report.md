# Bulow 2020 Analysis Notes

- Official paper-style figure entrypoints now live in `figure_2/plot_figure_2.py`, `figure_3/plot_figure_3.py`, `figure_4/plot_figure_4.py`, `figure_5/plot_figure_5.py`, `figure_6a/plot_figure_6a.py`, and `figure_6b/plot_figure_6b.py`.
- Figures 2-7 now stage their canonical Gibbs-energy and MIAC inputs under each figure's `input/` directory instead of reading from the repo root.
- Figure 5 includes a zero-height Debye-Huckel bar because the supplied digitized contribution tables only contain `hc`, `disp`, `assoc`, and `Born`.
- The diagnostic scans and non-paper auxiliary plots were kept under `diagnostics/` rather than deleted.
- The diagnostics directory still contains intermediate CSV and PNG artifacts from the earlier 2020 discrepancy work.
