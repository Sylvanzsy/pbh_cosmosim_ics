#!/usr/bin/env python3
"""
plot_mf_tables.py -- check plot for the default mass-function tables.

Draws the two tables in examples/mass_functions/ in the style of Fig. 1 of
Magaraggia et al. 2026 (orange dashed = Hasinger 2020, blue solid = this work),
and, if the upstream x=[...]/y=[...] files are given, overplots their raw points.

    plot_mf_tables.py [--raw-thiswork FA2.txt] [--raw-hasinger Hasinger.txt] [-o out.png]
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from import_mf_xy import read_xy

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default=os.path.join(here, "..", "examples", "mass_functions"))
ap.add_argument("--raw-thiswork")
ap.add_argument("--raw-hasinger")
ap.add_argument("-o", "--out", default="mf_tables_check.png")
a = ap.parse_args()

BLUE, ORANGE, INK, MUTED = "#0C7BDC", "#E66100", "#1f2430", "#6b7280"
tw = np.loadtxt(os.path.join(a.dir, "mf_magaraggia2026_thiswork.txt"))
ha = np.loadtxt(os.path.join(a.dir, "mf_hasinger2020.txt"))
raw_tw = read_xy(a.raw_thiswork) if a.raw_thiswork else None
raw_ha = read_xy(a.raw_hasinger) if a.raw_hasinger else None

plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"wspace": 0.22})
for ax, (lo, hi) in zip(axs, [(1e-6, 3e10), (1e3, 3e10)]):
    if raw_ha is not None:
        inside = raw_ha[0] <= 1e8
        ax.loglog(raw_ha[0][inside], raw_ha[1][inside], "o", ms=4.5, mfc="none", mec=ORANGE, mew=1.1,
                  label="Hasinger.txt, source points")
        ax.loglog(raw_ha[0][~inside], raw_ha[1][~inside], "x", ms=6, color=MUTED, mew=1.3,
                  label=r"Hasinger.txt points above $10^8\,M_\odot$ (not in table)")
    if raw_tw is not None:
        ax.loglog(raw_tw[0], raw_tw[1], ".", ms=4.5, color=BLUE, alpha=0.55, label="FA2.txt, source points")
    ax.loglog(ha[:, 0], ha[:, 1], "--", color=ORANGE, lw=2.0, label="table: Hasinger 2020 (orange dashed in Fig. 1)")
    ax.loglog(tw[:, 0], tw[:, 1], "-", color=BLUE, lw=2.0, label="table: This work (blue solid in Fig. 1)")
    ax.axvspan(1e8, hi, color="#9ca3af", alpha=0.18, lw=0)
    ax.set_xlim(lo, hi)
    ax.set_xlabel(r"$M_{\rm PBH}\ [M_\odot]$")
    ax.grid(True, ls="--", lw=0.5, alpha=0.6)
axs[0].set_ylabel(r"$d f_{\rm PBH}/d\ln M$")
axs[0].set_ylim(1e-8, 1.0)
axs[0].set_title("Full range (tables end at $10^8\\,M_\\odot$)", color=INK)
axs[1].set_ylim(1e-8, 3e-4)
axs[1].set_title(r"Zoom: $M \geq 10^3\,M_\odot$", color=INK)
axs[0].text(1.6e8, 0.25, "outside\ntables", fontsize=9, color=MUTED, va="top")
axs[0].legend(frameon=False, fontsize=8.8, loc="lower left")
fig.savefig(a.out, dpi=140, bbox_inches="tight")
print("wrote", a.out)
