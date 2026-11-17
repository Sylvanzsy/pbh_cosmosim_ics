#!/usr/bin/env python3
"""
import_mf_xy.py -- convert a PBH mass function stored as two Python-list lines
(x = [...] masses in Msun, y = [...] df_PBH/dlnM), the format used by
https://github.com/albertomagaraggia/PBH_CosmicArchitects, into the two-column
table read by pbh_ic (MASSFUNC_TABLE, mf_type='dfdlnM'):

    # M [Msun]    df_PBH/dlnM

Usage:
    import_mf_xy.py IN.txt OUT.txt --label "..." [--source "..."]
    import_mf_xy.py --summary TABLE.txt        # print integrals of a table

By default the values are written as given. With --figure-curve the script
instead reproduces the curve drawn in Fig. 1 of Magaraggia et al. 2026: a
log-log smoothing spline through ALL points of the file
(scipy UnivariateSpline(log10 M, log10 psi, s=0.08), as in their notebook),
evaluated on 141 points from 1e-6 to 1e8 Msun (0.1 dex). Nothing is renormalised.
"""
import argparse, ast
import numpy as np


def read_xy(path):
    x = y = None
    for line in open(path).read().splitlines():
        s = line.lstrip()
        if s.startswith("x"):
            x = ast.literal_eval(s.split("=", 1)[1].strip().rstrip("."))
        elif s.startswith("y"):
            y = ast.literal_eval(s.split("=", 1)[1].strip().rstrip("."))
    if x is None or y is None or len(x) != len(y):
        raise SystemExit(f"{path}: could not read matching x=[...] and y=[...] lines")
    x, y = np.array(x, float), np.array(y, float)
    o = np.argsort(x)
    return x[o], y[o]


def f_between(M, psi, a, b, n=20000):
    """int psi dlnM over [a,b], with ln(psi) linear in ln(M) between table points."""
    a, b = max(a, M[0]), min(b, M[-1])
    if b <= a:
        return 0.0
    g = np.linspace(np.log(a), np.log(b), n)
    p = np.exp(np.interp(g, np.log(M), np.log(psi)))
    return float(np.sum(0.5 * (p[1:] + p[:-1]) * np.diff(g)))


def mean_mass(M, psi, a, b, n=20000):
    """Number-weighted mean mass in [a,b]: f / int (psi/M) dlnM."""
    a, b = max(a, M[0]), min(b, M[-1])
    g = np.linspace(np.log(a), np.log(b), n)
    p = np.exp(np.interp(g, np.log(M), np.log(psi)))
    w = p / np.exp(g)
    return f_between(M, psi, a, b, n) / float(np.sum(0.5 * (w[1:] + w[:-1]) * np.diff(g)))


def summary(M, psi):
    print(f"  points: {len(M)}   mass range: {M[0]:.3e} - {M[-1]:.3e} Msun")
    print(f"  total f (whole table)       = {f_between(M, psi, M[0], M[-1]):.4g}")
    for a, b in [(1e-6, 1e8), (1e3, 1e8), (1e3, 1e5), (1e5, 1e7), (1e5, 1e8), (1e6, 1e8)]:
        f = f_between(M, psi, a, b)
        if f > 0:
            print(f"  f in [{a:.0e}, {b:.0e}] Msun = {f:.4g}   <M>_number = {mean_mass(M, psi, a, b):.4g} Msun")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("inp", nargs="?")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--label", default="")
    ap.add_argument("--source", default="")
    ap.add_argument("--summary")
    ap.add_argument("--figure-curve", action="store_true")
    ap.add_argument("--smooth", type=float, default=0.08)
    a = ap.parse_args()
    if a.summary:
        d = np.loadtxt(a.summary, comments="#")
        summary(d[:, 0], d[:, 1])
        raise SystemExit
    M, psi = read_xy(a.inp)
    keep = (M > 0) & (psi > 0)
    M, psi = M[keep], psi[keep]
    how = "Values are copied unchanged from the source file (no smoothing, no renormalisation)."
    if a.figure_curve:
        from scipy.interpolate import UnivariateSpline
        spl = UnivariateSpline(np.log10(M), np.log10(psi), s=a.smooth)
        M = np.logspace(-6, 8, 141)
        psi = 10.0 ** spl(np.log10(M))
        how = (f"Values are the curve drawn in Fig. 1 of Magaraggia et al. 2026: log-log smoothing spline "
               f"(s={a.smooth:g}) through all source points, evaluated from 1e-6 to 1e8 Msun in 0.1 dex steps. "
               "Not renormalised.")
    hdr = [a.label, a.source,
           "Convention: second column is psi(M) = d f_PBH / d ln M (mass-weighted), M in Msun.",
           how,
           "For pbh_ic: MASSFUNC_TABLE with mf_type = 'dfdlnM'.",
           "M [Msun]    df_PBH/dlnM"]
    np.savetxt(a.out, np.c_[M, psi], fmt="%.10e", header="\n".join(h for h in hdr if h))
    summary(M, psi)
