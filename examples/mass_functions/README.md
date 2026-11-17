# Default extended PBH mass functions

Two-column tables, `M [Msun]` and `psi(M) = d f_PBH / d ln M` (mass-weighted).
They are the input for the tabulated mass-function mode of `pbh_ic`
(`MASSFUNC_TABLE`, `mf_type = 'dfdlnM'`). That mode is not in the code yet;
these files only provide the data.

| File | Mass function | Table range | Total f in table |
|---|---|---|---|
| `mf_magaraggia2026_FA2.txt` | Magaraggia et al. 2026 (arXiv:2608.26247), "this work": thermal history with lepton-flavour asymmetry | 1e-6 to 1e8 Msun, 141 points | 0.339 |
| `mf_hasinger2020.txt` | Hasinger 2020 thermal-history spectrum, as adopted by Cappelluti, Hasinger & Natarajan 2022 (doi:10.3847/1538-4357/ac332d) | 1.4e-8 to 6.7e9 Msun, 66 points | 0.872 |

Both tables are copied unchanged from
https://github.com/albertomagaraggia/PBH_CosmicArchitects (commit 9734cb3, files
`FA2.txt` and `Hasinger.txt`, MIT License, see `LICENSE_PBH_CosmicArchitects`).
The ApJ paper shows its spectrum only as a figure, so `mf_hasinger2020.txt` is the
Hasinger (2020) reference curve distributed with Magaraggia et al., not data taken
from the ApJ paper. Cite the original papers when using them.

Mass fraction of dark matter in the ranges relevant for simulations
(`ln psi` interpolated linearly in `ln M`):

| Range [Msun] | f (Magaraggia FA2) | f (Hasinger 2020) | number-weighted mean mass |
|---|---|---|---|
| 1e3 to 1e5 | 2.15e-4 | 2.13e-4 | 2.7e3 |
| 1e5 to 1e7 | 3.94e-5 | 3.78e-5 | 3.8e5 |
| 1e6 to 1e8 | 1.94e-5 | 1.90e-5 | 2.7e6 |
| 1e3 to 1e8 | 2.58e-4 | 2.53e-4 | 3.3e3 |

Above 1e3 Msun the two tables agree to within 5%. They differ at lower masses.

`tools/import_mf_xy.py` regenerates the tables from the upstream files and prints
these integrals (`--summary TABLE.txt`).
