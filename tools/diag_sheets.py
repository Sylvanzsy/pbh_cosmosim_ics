#!/usr/bin/env python3
"""Diagnose periodic-image sheet artefacts in PHANTOM's PBH perturbation.

Reads the PHANTOM output IC (post-perturbation) together with the original
MUSIC IC (pre-perturbation), matches particles by ID, computes the applied
displacement dx = x_out - x_in (minimum image), deposits div(dx) on a CIC
grid, and reports where div has sharp negative spikes. Those are compared
against:
  (a) x_PBH +/- L/2 and y/z_PBH +/- L/2 for every PBH (read from the
      <icbase>_pbh_list.txt catalogue) and axis - the minimum-image
      discontinuity predicted in the issue report;
  (b) multiples of L/nx, if --nx is given (PBH sampling-cell boundaries);
  (c) the rpert truncation radius, if --rpert is given (a sphere, so it
      would NOT show up as a planar feature - listed only to rule it out).

Also saves a 2D slice plot of div(dx) through the most massive PBH, with
the predicted +/-L/2 planes overlaid.

Usage:
    diag_sheets.py <output_snapshot> <input_ic> <pbh_catalogue>
        [--grid 128] [--nx NX] [--rpert RPERT] [--out sheets.png]
"""
import argparse
import sys

import numpy as np

HEADER_FIELDS = [
    ('npart', '<i4', 6), ('massarr', '<f8', 6), ('a', '<f8', 1),
    ('redshift', '<f8', 1), ('flag_sfr', '<i4', 1), ('flag_feedback', '<i4', 1),
    ('nall', '<i4', 6), ('flag_cool', '<i4', 1), ('nfs', '<i4', 1),
    ('boxsize', '<f8', 1), ('Om0', '<f8', 1), ('Ol', '<f8', 1), ('h', '<f8', 1),
    ('unused', '<i4', 24),
]
HEADER_SIZE = 256


def read_record(f):
    head = f.read(4)
    if len(head) < 4:
        return None
    n1 = int(np.frombuffer(head, dtype='<i4', count=1)[0])
    data = f.read(n1)
    n2 = int(np.frombuffer(f.read(4), dtype='<i4', count=1)[0])
    if n1 != n2:
        raise ValueError(f'record marker mismatch: {n1} != {n2}')
    return data


def parse_header(buf):
    off = 0
    hdr = {}
    for name, dt, count in HEADER_FIELDS:
        itemsize = np.dtype(dt).itemsize
        arr = np.frombuffer(buf, dtype=dt, count=count, offset=off)
        off += itemsize * count
        hdr[name] = arr.copy() if count > 1 else arr[0]
    return hdr


def read_gadget1(path):
    with open(path, 'rb') as f:
        hdr = parse_header(read_record(f))
        npart = hdr['npart']
        N = int(npart.sum())
        pos = np.frombuffer(read_record(f), dtype='<f4').reshape(N, 3).astype(np.float64)
        vel = np.frombuffer(read_record(f), dtype='<f4').reshape(N, 3).astype(np.float64)
        ids = np.frombuffer(read_record(f), dtype='<i4').copy()
    return dict(hdr=hdr, pos=pos, vel=vel, ids=ids)


def load_catalogue(path):
    ids, pos, mass = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            tok = line.split()
            ids.append(int(tok[1]))
            pos.append([float(x) for x in tok[2:5]])
            mass.append(float(tok[8]))
    return (np.array(ids, dtype=np.int64), np.array(pos, dtype=np.float64),
            np.array(mass, dtype=np.float64))


def cic_deposit_vector(pos, vec, Ng, bs):
    """CIC-deposit a per-particle vector field onto an Ng^3 grid, returning
    the particle-count-weighted average field per cell (and the count)."""
    cell = bs / Ng
    g = pos / cell
    i0 = np.floor(g).astype(np.int64)
    frac = g - i0
    sumw = np.zeros((Ng, Ng, Ng), dtype=np.float64)
    sumv = np.zeros((Ng, Ng, Ng, 3), dtype=np.float64)
    for dx in (0, 1):
        wx = (1 - frac[:, 0]) if dx == 0 else frac[:, 0]
        ix = (i0[:, 0] + dx) % Ng
        for dy in (0, 1):
            wy = (1 - frac[:, 1]) if dy == 0 else frac[:, 1]
            iy = (i0[:, 1] + dy) % Ng
            for dz in (0, 1):
                wz = (1 - frac[:, 2]) if dz == 0 else frac[:, 2]
                iz = (i0[:, 2] + dz) % Ng
                w = wx * wy * wz
                flat = (ix * Ng + iy) * Ng + iz
                np.add.at(sumw.reshape(-1), flat, w)
                for c in range(3):
                    np.add.at(sumv.reshape(-1, 3)[:, c], flat, w * vec[:, c])
    psi = np.divide(sumv, sumw[..., None], out=np.zeros_like(sumv),
                     where=sumw[..., None] > 0)
    return psi, sumw


def periodic_divergence(psi, cell):
    dpx = (np.roll(psi[..., 0], -1, axis=0) - np.roll(psi[..., 0], 1, axis=0)) / (2 * cell)
    dpy = (np.roll(psi[..., 1], -1, axis=1) - np.roll(psi[..., 1], 1, axis=1)) / (2 * cell)
    dpz = (np.roll(psi[..., 2], -1, axis=2) - np.roll(psi[..., 2], 1, axis=2)) / (2 * cell)
    return dpx + dpy + dpz


def near_pbh_mask(Ng, cell, bs, pbh_pos, exclude_cells):
    """True for grid cells within `exclude_cells` grid spacings (periodic) of
    any PBH - the genuine, physically-expected near-field convergence around
    a point source, which must be masked out before profiling for the much
    weaker, box-spanning +/-L/2 sheet: otherwise the near-field spike (huge
    but localised to a few cells) dominates a plane average and hides it."""
    idx = np.arange(Ng)
    coord = idx * cell
    mask = np.zeros((Ng, Ng, Ng), dtype=bool)
    r_excl = exclude_cells * cell
    for p in pbh_pos:
        d = np.abs(coord[:, None] - p[None, :])  # (Ng,3): d[:,c] vs PBH component c
        d = np.minimum(d, bs - d)
        dx2 = (d[:, 0] ** 2).reshape(Ng, 1, 1)
        dy2 = (d[:, 1] ** 2).reshape(1, Ng, 1)
        dz2 = (d[:, 2] ** 2).reshape(1, 1, Ng)
        sep2 = dx2 + dy2 + dz2
        mask |= sep2 < r_excl ** 2
    return mask


def find_spike_planes(div, cell, axis, nsig=5.0, mask=None):
    """Plane-averaged divergence along `axis` (cells under `mask` excluded,
    e.g. the PBH near field); flag indices with a sharp spike - positive
    (underdense, div(dx)>0 i.e. particles pulled apart) or negative
    (overdense) - relative to the typical (robust) scatter of the
    remaining, far-field cells. div(dx)=-delta in the Zel'dovich
    approximation, so the sign that shows up depends on which side of the
    discontinuity dominates the plane average; both are reported."""
    other = tuple(a for a in range(3) if a != axis)
    d = div.copy()
    if mask is not None:
        d[mask] = np.nan
    prof = np.nanmean(d, axis=other)
    med = np.nanmedian(prof)
    mad = np.nanmedian(np.abs(prof - med)) * 1.4826 + 1e-300
    z = (prof - med) / mad
    idx = np.where(np.abs(z) > nsig)[0]
    order = np.argsort(-np.abs(z[idx]))
    idx = idx[order]
    return idx, idx * cell, prof, z


def nearest_match(value, candidates, bs):
    if len(candidates) == 0:
        return None, None
    d = np.abs(candidates - value)
    d = np.minimum(d, bs - d)
    k = np.argmin(d)
    return candidates[k], d[k]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('output_snapshot')
    ap.add_argument('input_ic')
    ap.add_argument('pbh_catalogue')
    ap.add_argument('--grid', type=int, default=128)
    ap.add_argument('--nx', type=int, default=None,
                     help='PBH sampling-grid size per side (optional cross-check)')
    ap.add_argument('--rpert', type=float, default=None,
                     help='rpert truncation radius [kpc/h] (optional cross-check)')
    ap.add_argument('--out', default='sheets.png')
    ap.add_argument('--nsig', type=float, default=5.0,
                     help='MAD-sigma threshold for flagging a plane spike')
    ap.add_argument('--exclude-cells', type=float, default=4.0,
                     help='mask out cells within this many grid spacings of any PBH '
                          'before profiling (excludes the genuine near-field '
                          'convergence spike so it does not hide the +/-L/2 sheet)')
    args = ap.parse_args()

    out = read_gadget1(args.output_snapshot)
    inp = read_gadget1(args.input_ic)
    bs = float(out['hdr']['boxsize'])
    Ng = args.grid
    cell = bs / Ng

    id_to_row_in = {int(i): r for r, i in enumerate(inp['ids'])}
    common_ids = np.array([i for i in out['ids'] if int(i) in id_to_row_in], dtype=np.int64)
    rows_out = {int(i): r for r, i in enumerate(out['ids'])}
    out_idx = np.array([rows_out[int(i)] for i in common_ids])
    in_idx = np.array([id_to_row_in[int(i)] for i in common_ids])
    print(f'Matched {len(common_ids)} particles by ID '
          f'(output has {len(out["ids"])}, input has {len(inp["ids"])})')

    pos_out = out['pos'][out_idx]
    pos_in = inp['pos'][in_idx]
    dx = pos_out - pos_in
    dx -= bs * np.round(dx / bs)  # minimum image
    print('max |dx| component [kpc/h] =', np.max(np.abs(dx)))

    psi, weight = cic_deposit_vector(pos_out, dx, Ng, bs)
    div = periodic_divergence(psi, cell)
    print(f'div(dx) on {Ng}^3 grid: min={div.min():.4e} max={div.max():.4e} '
          f'std={div.std():.4e}')

    pbh_ids, pbh_pos, pbh_mass = load_catalogue(args.pbh_catalogue)
    if len(pbh_mass) == 0:
        print('No PBHs in catalogue - nothing to compare against.')
        sys.exit(0)
    i_heaviest = int(np.argmax(pbh_mass))
    print(f'Heaviest PBH: id={pbh_ids[i_heaviest]} mass={pbh_mass[i_heaviest]:.4e} Msun '
          f'pos={pbh_pos[i_heaviest]}')

    near_mask = near_pbh_mask(Ng, cell, bs, pbh_pos, args.exclude_cells)
    print(f'Masking {near_mask.sum()}/{near_mask.size} cells within '
          f'{args.exclude_cells:.1f} grid spacings of a PBH before profiling.')

    axis_names = ['x', 'y', 'z']
    detected = {}
    for axis in range(3):
        idx, coords, prof, z = find_spike_planes(div, cell, axis, nsig=args.nsig, mask=near_mask)
        detected[axis] = coords
        print(f'\n--- axis {axis_names[axis]}: {len(coords)} spike plane(s) detected ---')
        for i, c in zip(idx, coords):
            sign = 'underdense (div>0)' if z[i] > 0 else 'overdense (div<0)'
            print(f'  plane at {axis_names[axis]} = {c:9.3f} kpc/h  '
                  f'(z={z[i]:+.2f}, {sign})')

    print('\n=== cross-check: detected planes vs x_PBH +/- L/2 (predicted by the '
          'minimum-image diagnosis) ===')
    any_match = False
    for p in range(len(pbh_mass)):
        for axis in range(3):
            for sign in (+1, -1):
                pred = (pbh_pos[p, axis] + sign * bs / 2.0) % bs
                match, d = nearest_match(pred, detected[axis], bs)
                tag = 'MATCH' if (match is not None and d < 2 * cell) else 'no match'
                if tag == 'MATCH':
                    any_match = True
                print(f'  PBH#{p} (M={pbh_mass[p]:.2e} Msun) axis={axis_names[axis]} '
                      f'predicted={pred:9.3f}  nearest detected={match}  '
                      f'diff={d if d is not None else float("nan"):.3f}  [{tag}]')
    print(f'\nAny detected plane matches a predicted x_PBH +/- L/2 plane: {any_match}')

    if args.nx is not None:
        cell_pbh = bs / args.nx
        grid_planes = np.arange(args.nx + 1) * cell_pbh
        print(f'\n=== cross-check: detected planes vs PBH sampling-cell boundaries '
              f'(multiples of L/{args.nx} = {cell_pbh:.3f} kpc/h) ===')
        for axis in range(3):
            for c in detected[axis]:
                match, d = nearest_match(c, grid_planes, bs)
                print(f'  axis={axis_names[axis]} detected={c:9.3f}  '
                      f'nearest sampling-cell boundary={match:.3f}  diff={d:.3f}')

    if args.rpert is not None:
        print(f'\nrpert={args.rpert:.3f} kpc/h: a truncation RADIUS produces a spherical '
              f'(not planar) feature around each PBH, so it would not explain a planar '
              f'spike at fixed x or y - listed for elimination only.')

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # A single-cell slice is dominated by per-cell shot noise - the sheet
        # only becomes visually obvious the same way it was found
        # numerically, by averaging along the third axis (like the z=47
        # column projection that first showed it). Project over z here.
        proj_raw = np.nanmean(div, axis=2).T
        proj_far = np.nanmean(np.where(near_mask, np.nan, div), axis=2).T
        far_scale = np.nanpercentile(np.abs(proj_far), 99.5) + 1e-300

        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        panels = [
            (axes[0], proj_raw, np.nanmax(np.abs(proj_raw)),
             'full range, z-projected (near-field PBH column dominates)'),
            (axes[1], proj_far, far_scale,
             'near field masked (reveals the +/-L/2 sheet)'),
        ]
        for ax, sl, vmax, title in panels:
            im = ax.imshow(sl, origin='lower', extent=[0, bs, 0, bs],
                            cmap='RdBu_r', vmin=-vmax, vmax=vmax)
            for sign in (+1, -1):
                ax.axvline((pbh_pos[i_heaviest, 0] + sign * bs / 2.0) % bs,
                           color='k', ls='--', lw=1, alpha=0.7)
                ax.axhline((pbh_pos[i_heaviest, 1] + sign * bs / 2.0) % bs,
                           color='k', ls='--', lw=1, alpha=0.7)
            ax.plot(pbh_pos[i_heaviest, 0], pbh_pos[i_heaviest, 1], 'k*', ms=14)
            ax.set_xlabel('x [kpc/h]')
            ax.set_ylabel('y [kpc/h]')
            ax.set_title(title, fontsize=9)
            fig.colorbar(im, ax=ax, label='div(dx)')
        fig.suptitle('div(dx), projected over z, through the heaviest PBH '
                     '(dashed: predicted x_PBH,y_PBH +/- L/2)')
        fig.tight_layout()
        fig.savefig(args.out, dpi=150)
        print(f'\nSaved slice plot to {args.out}')
    except ImportError:
        print('\nmatplotlib not available - skipping plot.')


if __name__ == '__main__':
    main()
