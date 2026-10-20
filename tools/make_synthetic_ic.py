#!/usr/bin/env python3
"""Build a tiny synthetic single-file Gadget format-1 IC (DM-only).

Matches the exact header layout read by pbh_ic.F90:

    read (1) npart, massarr, a, redshift, flag_sfr, flag_feedback, nall,
             flag_cool, nfs, boxsize, Om0, Ol, h, unused

followed by pos, vel, id. massarr(1) is set non-zero (fixed-mass DM, no
variable-mass block) and npart(0)=0 (no gas block), so the body is just
three records: pos, vel, id - the minimal case pbh_ic.F90 supports.

Used for a cheap end-to-end test of the pbh_ic.F90 changes without running
the real production-size IC. See tools/check_snapshot.py for the matching
read-back checker.
"""
import argparse

import numpy as np


def write_record(f, arr):
    data = np.ascontiguousarray(arr).tobytes()
    n = len(data)
    marker = np.array([n], dtype='<i4').tobytes()
    f.write(marker)
    f.write(data)
    f.write(marker)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('out', help='output <icbase>_cdm.dat path')
    ap.add_argument('--n1d', type=int, default=8, help='DM particles per side (N=n1d**3)')
    ap.add_argument('--boxsize', type=float, default=1000.0, help='box size [kpc/h]')
    ap.add_argument('--mdm', type=float, default=1.0e-3, help='DM particle mass [1e10 Msun/h]')
    ap.add_argument('--h', type=float, default=0.6774)
    ap.add_argument('--om0', type=float, default=0.3089)
    ap.add_argument('--ol', type=float, default=0.6911)
    ap.add_argument('--redshift', type=float, default=99.0)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    Ndm = args.n1d ** 3
    a = 1.0 / (1.0 + args.redshift)

    npart = np.zeros(6, dtype='<i4')
    npart[1] = Ndm
    nall = npart.copy()
    massarr = np.zeros(6, dtype='<f8')
    massarr[1] = args.mdm

    header = b''.join([
        npart.tobytes(),
        massarr.tobytes(),
        np.array(a, dtype='<f8').tobytes(),
        np.array(args.redshift, dtype='<f8').tobytes(),
        np.array(0, dtype='<i4').tobytes(),  # flag_sfr
        np.array(0, dtype='<i4').tobytes(),  # flag_feedback
        nall.tobytes(),
        np.array(0, dtype='<i4').tobytes(),  # flag_cool
        np.array(1, dtype='<i4').tobytes(),  # nfs (NumFiles)
        np.array(args.boxsize, dtype='<f8').tobytes(),
        np.array(args.om0, dtype='<f8').tobytes(),
        np.array(args.ol, dtype='<f8').tobytes(),
        np.array(args.h, dtype='<f8').tobytes(),
        np.zeros(24, dtype='<i4').tobytes(),  # unused
    ])
    assert len(header) == 256, len(header)

    # Fortran pos(1:3,1:N) is column-major -> (x,y,z) contiguous per particle,
    # i.e. the same byte layout as a C-order numpy (N,3) array flattened.
    pos = (rng.random((Ndm, 3)).astype('<f4')) * np.float32(args.boxsize)
    vel = rng.normal(0.0, 50.0, (Ndm, 3)).astype('<f4')
    ids = np.arange(Ndm, dtype='<i4')

    with open(args.out, 'wb') as f:
        write_record(f, np.frombuffer(header, dtype='<u1'))
        write_record(f, pos.reshape(-1))
        write_record(f, vel.reshape(-1))
        write_record(f, ids)

    print(f'Wrote {args.out}: Ndm={Ndm}, boxsize={args.boxsize} kpc/h, '
          f'mdm={args.mdm:.6e} (1e10 Msun/h), h={args.h}')


if __name__ == '__main__':
    main()
