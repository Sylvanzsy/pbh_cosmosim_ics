# Default extended PBH mass functions

Two-column tables, `M [Msun]` and `psi(M) = d f_PBH / d ln M` (mass-weighted).
They are the input for the tabulated mass-function mode of `pbh_ic`
(`MASSFUNC_TABLE`, `mf_type = 'dfdlnM'`). That mode is not in the code yet;
these files only provide the data.

Both tables are the two curves of Fig. 1 of Magaraggia et al. 2026
(arXiv:2608.26247), on the range plotted there, 1e-6 to 1e8 Msun (141 points,
0.1 dex steps):

| File | Curve in Fig. 1 | Total f in table |
|---|---|---|
| `mf_magaraggia2026_thiswork.txt` | blue solid, "This work": thermal history with lepton-flavour asymmetry | 0.338 |
| `mf_hasinger2020.txt` | orange dashed, "Hasinger 2020": the spectrum used by Cappelluti, Hasinger & Natarajan 2022 (doi:10.3847/1538-4357/ac332d) | 0.730 |

How they were made: the source points are the files `FA2.txt` and `Hasinger.txt`
of https://github.com/albertomagaraggia/PBH_CosmicArchitects (commit 9734cb3,
MIT License, see `LICENSE_PBH_CosmicArchitects`). As in that paper's notebook, a
log-log smoothing spline (`UnivariateSpline`, s = 0.08) is put through all source
points and evaluated on the grid above. The tables follow the source points to
about 6% rms (13% at most) above 1e3 Msun. `Hasinger.txt` also holds seven points
between 1e8 and 7e9 Msun; these are outside the published curve and are not in
the table. Nothing is renormalised. Cite the original papers when using them.

`mf_tables_check.png` shows both tables with the source points.

Mass fraction of dark matter in the ranges relevant for simulations
(`ln psi` interpolated linearly in `ln M`):

| Range [Msun] | f (This work) | f (Hasinger 2020) | number-weighted mean mass |
|---|---|---|---|
| 1e3 to 1e5 | 2.09e-4 | 2.10e-4 | 2.8e3 |
| 1e5 to 1e7 | 3.94e-5 | 3.80e-5 | 3.8e5 |
| 1e6 to 1e8 | 1.90e-5 | 1.93e-5 | 2.7e6 |
| 1e3 to 1e8 | 2.51e-4 | 2.51e-4 | 3.4e3 |

Above 1e3 Msun the two curves are nearly the same. They differ at lower masses.

Tools: `tools/import_mf_xy.py --figure-curve` regenerates the tables from the
upstream files and prints these integrals (`--summary TABLE.txt`);
`tools/plot_mf_tables.py` redraws the check plot.
