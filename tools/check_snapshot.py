#!/usr/bin/env python3
"""Read-back checker for pbh_ic.F90 output snapshots.

Verifies, using only header-derived block sizes (no magic offsets):
  1. The type-pbh_type (default 3) particles in the output snapshot match
     the <icbase>_pbh_list.txt catalogue row-by-row, matched by ID.
  2. Total mass is conserved: DM mass (output, post-reduction) + PBH mass
     (output) equals the DM mass in the INPUT ic (pre-PBH), to high
     relative precision.

Usage:
    check_snapshot.py <output_snapshot.dat> <pbh_catalogue.txt> <input_ic.dat>

Exits nonzero on any check failure.

Caveat on the total-mass tolerance: when the snapshot's PBH masses live in
the real*4 variable-mass block (MASSFUNC_LOGNORMAL path), the achievable
precision is bounded by float32 storage (~1e-7 relative), not the ~1e-12
that holds when PBH mass is monochromatic and stored as a real*8 header
value (massarr(pbh_type)). Pass --tol-total-mass to relax the default.
"""
import argparse
import sys

import numpy as np

HEADER_FIELDS = [
    ('npart', '<i4', 6),
    ('massarr', '<f8', 6),
    ('a', '<f8', 1),
    ('redshift', '<f8', 1),
    ('flag_sfr', '<i4', 1),
    ('flag_feedback', '<i4', 1),
    ('nall', '<i4', 6),
    ('flag_cool', '<i4', 1),
    ('nfs', '<i4', 1),
    ('boxsize', '<f8', 1),
    ('Om0', '<f8', 1),
    ('Ol', '<f8', 1),
    ('h', '<f8', 1),
    ('unused', '<i4', 24),
]
HEADER_SIZE = 256


def read_record(f):
    """Read one gfortran unformatted-sequential record (4-byte markers)."""
    head = f.read(4)
    if len(head) < 4:
        return None
    n1 = int(np.frombuffer(head, dtype='<i4', count=1)[0])
    data = f.read(n1)
    if len(data) != n1:
        raise ValueError(f'short read: expected {n1} bytes, got {len(data)}')
    n2 = int(np.frombuffer(f.read(4), dtype='<i4', count=1)[0])
    if n1 != n2:
        raise ValueError(f'record marker mismatch: leading={n1} trailing={n2} '
                          '(not a 4-byte-marker gfortran unformatted file?)')
    return data


def parse_header(buf):
    if len(buf) != HEADER_SIZE:
        raise ValueError(f'header record is {len(buf)} bytes, expected {HEADER_SIZE}')
    off = 0
    hdr = {}
    for name, dt, count in HEADER_FIELDS:
        itemsize = np.dtype(dt).itemsize
        arr = np.frombuffer(buf, dtype=dt, count=count, offset=off)
        off += itemsize * count
        hdr[name] = arr.copy() if count > 1 else arr[0]
    assert off == HEADER_SIZE
    return hdr


def mass_block_offsets(npart, massarr):
    """Offset into the variable-mass block for each type with massarr==0 and npart>0,
    in type order - matches Fortran's 'Nm = Nm + nall(i)' accumulation."""
    offs = {}
    running = 0
    for t in range(6):
        if massarr[t] == 0.0 and npart[t] > 0:
            offs[t] = running
            running += int(npart[t])
    return offs, running


def read_gadget1(path):
    with open(path, 'rb') as f:
        hdr = parse_header(read_record(f))
        npart = hdr['npart']
        N = int(npart.sum())
        pos = np.frombuffer(read_record(f), dtype='<f4').reshape(N, 3)
        vel = np.frombuffer(read_record(f), dtype='<f4').reshape(N, 3)
        ids = np.frombuffer(read_record(f), dtype='<i4').copy()
        mass_offs, Nm = mass_block_offsets(npart, hdr['massarr'])
        if Nm > 0:
            lmass = np.frombuffer(read_record(f), dtype='<f4').copy()
        else:
            lmass = np.zeros(0, dtype='<f4')
        Ngas = int(npart[0])
        if Ngas > 0:
            iegas = np.frombuffer(read_record(f), dtype='<f4').copy()
        else:
            iegas = np.zeros(0, dtype='<f4')
    return dict(hdr=hdr, pos=pos, vel=vel, ids=ids, lmass=lmass,
                mass_offs=mass_offs, iegas=iegas)


def type_slice(npart, t):
    start = int(npart[:t].sum())
    return slice(start, start + int(npart[t]))


def masses_for_type(snap, t):
    npart = snap['hdr']['npart']
    massarr = snap['hdr']['massarr']
    n_t = int(npart[t])
    if massarr[t] != 0.0:
        return np.full(n_t, massarr[t], dtype=np.float64)
    off = snap['mass_offs'][t]
    return snap['lmass'][off:off + n_t].astype(np.float64)


def load_catalogue(path):
    ids, pos, vel, mass = [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            tok = line.split()
            ids.append(int(tok[1]))
            pos.append([float(x) for x in tok[2:5]])
            vel.append([float(x) for x in tok[5:8]])
            mass.append(float(tok[8]))
    return (np.array(ids, dtype=np.int64), np.array(pos, dtype=np.float64),
            np.array(vel, dtype=np.float64), np.array(mass, dtype=np.float64))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('output_snapshot')
    ap.add_argument('pbh_catalogue')
    ap.add_argument('input_ic')
    ap.add_argument('--pbh-type', type=int, default=3)
    ap.add_argument('--tol-mass', type=float, default=1e-5,
                     help='max relative mass error, snapshot vs catalogue (float32-bounded)')
    ap.add_argument('--tol-pos', type=float, default=1e-2,
                     help='max |position difference| [kpc/h], snapshot vs catalogue')
    ap.add_argument('--tol-total-mass', type=float, default=1e-9,
                     help='max relative mismatch in DM+PBH vs input DM mass')
    args = ap.parse_args()

    out = read_gadget1(args.output_snapshot)
    inp = read_gadget1(args.input_ic)

    npart = out['hdr']['npart']
    massarr = out['hdr']['massarr']
    h = float(out['hdr']['h'])
    pt = args.pbh_type
    n_pbh = int(npart[pt])
    sl = type_slice(npart, pt)

    pos_pbh = out['pos'][sl].astype(np.float64)
    vel_pbh = out['vel'][sl].astype(np.float64)
    id_pbh = out['ids'][sl]
    mass_pbh_msun = masses_for_type(out, pt) * 1.0e10 / h

    cat_ids, cat_pos, cat_vel, cat_mass = load_catalogue(args.pbh_catalogue)

    ok = True

    if len(cat_ids) != n_pbh:
        print(f'FAIL: catalogue has {len(cat_ids)} rows, snapshot has '
              f'{n_pbh} type-{pt} particles')
        ok = False

    order = np.argsort(id_pbh)
    cat_order = np.argsort(cat_ids)
    if not np.array_equal(id_pbh[order], cat_ids[cat_order]):
        print('FAIL: ID sets differ between snapshot and catalogue')
        ok = False
    else:
        pos_s, vel_s, mass_s = pos_pbh[order], vel_pbh[order], mass_pbh_msun[order]
        pos_c, vel_c, mass_c = cat_pos[cat_order], cat_vel[cat_order], cat_mass[cat_order]

        dpos = float(np.max(np.abs(pos_s - pos_c))) if n_pbh else 0.0
        dvel = float(np.max(np.abs(vel_s - vel_c))) if n_pbh else 0.0
        rel_mass = float(np.max(np.abs(mass_s - mass_c) / np.maximum(mass_c, 1e-300))) if n_pbh else 0.0

        print(f'n_pbh                          = {n_pbh}')
        print(f'max |pos diff|   [kpc/h]       = {dpos:.3e}  (tol {args.tol_pos:.1e})')
        print(f'max |vel diff|   [u]           = {dvel:.3e}')
        print(f'max relative mass error        = {rel_mass:.3e}  (tol {args.tol_mass:.1e})')

        if dpos > args.tol_pos:
            print('FAIL: position mismatch exceeds tolerance')
            ok = False
        if rel_mass > args.tol_mass:
            print('FAIL: mass mismatch exceeds tolerance')
            ok = False

    dm_in = float(inp['hdr']['npart'][1]) * float(inp['hdr']['massarr'][1]) * 1.0e10 / float(inp['hdr']['h'])
    dm_out = float(npart[1]) * float(massarr[1]) * 1.0e10 / h
    pbh_out = float(mass_pbh_msun.sum())
    rel_tot = abs((dm_out + pbh_out) - dm_in) / dm_in if dm_in != 0 else float('nan')

    print(f'DM mass (input)   [Msun]       = {dm_in:.8e}')
    print(f'DM mass (output)  [Msun]       = {dm_out:.8e}')
    print(f'PBH mass (output) [Msun]       = {pbh_out:.8e}')
    print(f'relative mismatch (DM_out+PBH_out - DM_in)/DM_in = {rel_tot:.3e}  '
          f'(tol {args.tol_total_mass:.1e})')
    if not (rel_tot <= args.tol_total_mass):
        print('FAIL: total mass not conserved within tolerance')
        ok = False

    print('CHECK PASSED' if ok else 'CHECK FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
