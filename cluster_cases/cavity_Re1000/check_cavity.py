#!/usr/bin/env python3
"""Lid-driven cavity Re = 1000 scene check from a bin/mixnsolver output .vtu.

Usage:  python3 check_cavity.py <case.vtu>

Case: unit square cavity [0,1] x [0,1], lid (y = 1) moving at U_lid = (1,0,0),
Re = rho*U*L/mu = 1*1*1/0.001 = 1000, strictly 2D (one cell in z).

What is printed
  * grid summary (solver cells, cell-centroid range, max|u|)
  * the u_x velocity on the VERTICAL centreline (x = 0.5) at the 17 Ghia
    stations, against the Ghia, Ghia & Shin (1982) Re = 1000 table, plus the
    relative L2 error over the interior stations (the two wall rows are
    excluded -- a cell-centred probe cannot reproduce the boundary value).
    This is exactly the comparison the solver itself prints (ghia_compare in
    src/unstructured/mod_uns_output.f90) when it runs the MPI path; this script
    works on any output, serial or parallel.

VTU layout of this code (see vtk_write in src/unstructured/mod_uns_output.f90):
points = mesh nodes then one cell-centroid point per cell; cells = tetrahedra
from a face fan + the cell centroid, written as "n_pts  i1 i2 i3  centroid_id"
(5 integers per tet).  Cell data are replicated on every sub-tet, so the unique
centroid ids recover one value per solver cell.
"""
import re
import sys

import numpy as np

# --- Ghia, Ghia & Shin (1982), JCP 48:387-411, Table I (Re = 1000) ----------
YG = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172,
      0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]
U1000 = [1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719,
         0.05702, -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220,
         -0.20196, -0.18109, 0.00000]


def read_vtu(fn):
    txt = open(fn).read()

    def arr(name, cols=1, dtype=float):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % re.escape(name),
                      txt, re.S)
        if m is None:
            raise SystemExit('check_cavity: no "%s" array in %s' % (name, fn))
        a = np.array([dtype(x) for x in m.group(1).split()])
        return a.reshape(-1, cols) if cols > 1 else a

    pts = arr('Coordinates', 3)
    n_cells = int(re.search(r'NumberOfCells="(\d+)"', txt).group(1))
    con = arr('connectivity', 5, int)
    if con.shape[0] != n_cells:
        raise SystemExit('check_cavity: unexpected connectivity length')
    vel = arr('velocity', 3)

    # one value per solver cell (first sub-tet of each centroid group)
    cen_id = con[:, 4]
    first = np.concatenate(([True], cen_id[1:] != cen_id[:-1]))
    idx = np.nonzero(first)[0]
    base = cen_id.min()
    xyz = pts[cen_id[idx], :]                  # cell centroids
    return xyz, vel[idx, :], cen_id[idx] - base + 1


def probe_velocity(xyz, vel, xq, knn=8):
    """Inverse-square-distance probe over the `knn` nearest cell centroids.

    Replicates probe_velocity in src/unstructured/mod_uns_output.f90 (KNN = 8,
    w = 1/(d^2 + 1e-12)) so the numbers below match the solver's own
    ghia_compare table exactly.
    """
    d2 = np.sum((xyz - xq) ** 2, axis=1)
    idx = np.argsort(d2, kind='stable')[:knn]
    w = 1.0 / (d2[idx] + 1.0e-12)
    return (w[:, None] * vel[idx]).sum(axis=0) / w.sum()


def main(fn):
    xyz, vel, cid = read_vtu(fn)
    print('file            : %s' % fn)
    print('solver cells    : %d' % cid.max())
    print('centroid range  : x %.4f..%.4f  y %.4f..%.4f  z %.4f..%.4f'
          % (xyz[:, 0].min(), xyz[:, 0].max(), xyz[:, 1].min(),
             xyz[:, 1].max(), xyz[:, 2].min(), xyz[:, 2].max()))
    print('max|u| / U_lid  : %.6f    (cell-centred max; the lid value 1.0 is a\n'
          '                              FACE value, so the first cell below it\n'
          '                              reads ~0.91 on this mesh)'
          % np.hypot(vel[:, 0], vel[:, 1]).max())

    print('')
    print('--- Ghia et al. (1982) comparison, Re = 1000, vertical centreline ---')
    print('      y        u_num      u_ghia      error')
    l2num = l2den = 0.0
    for i, (yg, ug) in enumerate(zip(YG, U1000), start=1):
        un = probe_velocity(xyz, vel, np.array([0.5, yg, 0.5]))[0]
        err = un - ug
        print('%9.4f%12.5e%12.5e%12.2e' % (yg, un, ug, err))
        if 1 < i < len(YG):                    # skip the two wall rows
            l2num += err * err
            l2den += ug * ug
    l2 = np.sqrt(l2num / max(l2den, np.finfo(float).tiny))
    print('relative L2 error (vertical-centreline u_x): %.5e' % l2)
    print('(the MPI solver prints the same quantity as '
          '"relative L2 error (centerline u_x)" in logs/<tag>.log)')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: check_cavity.py <case.vtu>')
    main(sys.argv[1])
