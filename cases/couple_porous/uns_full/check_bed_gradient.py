#!/usr/bin/env python3
"""Verify the uns_full reference reproduces the analytic bed pressure drop.

Reads the solver VTU written by mod_uns_output.f90 (vtk_write), whose layout
is: points = the 6363 mesh nodes followed by one CENTROID point per cell;
cells = 12 tetrahedra per hex (face fan + centroid apex), grouped cell by
cell, with every CellData value replicated over the cell's sub-tets.  The
5th node of each tet is therefore the centroid point id = nnodes + c - 1,
which lets us recover the exact per-cell value (a simple nodal average would
smear it across the fluid/porous interface).

The case is strictly 1D (lateral faces are symmetry), so the expected answer
is: flat pressure over the fluid half x=0..100 mm, then a constant Darcy +
Forchheimer gradient over the bed half x=100..200 mm:

    dP/dx = (mu/eps/K) u + rho beta u^2
          = 1602.2 + 1419.5 = 3021.7 Pa/m  ->  3.0217 Pa/mm
    dP    = 302.1 Pa over L = 0.1 m

Usage:  python3 check_bed_gradient.py [unMesh.vtu]
"""
import sys
import re
import numpy as np

NNODES = 6363                 # mesh nodes, from the generator (101*21*3)
L_BED = 0.1                   # m
DP_ANALYTIC = 302.1           # Pa
GRAD_ANALYTIC = 3021.7        # Pa/m


def load(path):
    txt = open(path).read()

    def arr(name):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        return np.fromstring(m.group(1), sep=' ')

    pts = arr('Coordinates').reshape(-1, 3)
    conn = arr('connectivity').astype(np.int64)
    pres = arr('pressure')
    # 5th node of every tet is its cell's centroid point id
    cen = conn[4::5]
    starts = np.flatnonzero(np.concatenate(([True], cen[1:] != cen[:-1])))
    cid0 = cen[starts] - NNODES
    order = np.argsort(cid0)
    cid0, starts = cid0[order], starts[order]
    return pts[NNODES + cid0, 0], pres[starts], cid0


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'unMesh.vtu'
    xc, p, cid0 = load(path)

    fluid = xc[cid0 < 2000]
    bed = xc[cid0 >= 2000]
    print('%s : %d cells recovered' % (path, p.size))
    print('  zone 2 (fluid) x = %.4f..%.4f m' % (fluid.min(), fluid.max()))
    print('  zone 8 (bed)   x = %.4f..%.4f m' % (bed.min(), bed.max()))
    assert fluid.max() <= 0.1001, 'zone 2 is not the x<100 mm half!'
    assert bed.min() >= 0.0999, 'zone 8 is not the x>100 mm half!'

    def fit(x0, x1):
        m = (xc > x0) & (xc < x1)
        A = np.vstack([xc[m], np.ones(int(m.sum()))]).T
        return np.linalg.lstsq(A, p[m], rcond=None)[0]

    # fit windows exclude the first/last few planes, which carry the
    # mass-flow-inlet ramp and pressure-outlet "adjustment" artefacts
    # (see README.md); the physical gradient is exact in between.
    g_fluid, _ = fit(0.005, 0.095)
    g_bed, c_bed = fit(0.105, 0.195)
    # plane-averaged pressures at the mesh planes straddling the interface
    p_face = p[(xc > 0.098) & (xc < 0.100)].mean()   # last fluid plane
    p_bed0 = p[(xc > 0.100) & (xc < 0.102)].mean()   # first bed plane
    p_out = p[xc >= 0.198].mean()                    # last bed plane
    print()
    print('  fluid half gradient : %9.1f Pa/m  (%.4f Pa/mm, expect ~0)'
          % (g_fluid, g_fluid / 1000.0))
    print('  bed   half gradient : %9.1f Pa/m  (%.4f Pa/mm, expect %.4f)'
          % (g_bed, g_bed / 1000.0, GRAD_ANALYTIC / 1000.0))

    dp_bed = -g_bed * L_BED
    print()
    print('  p (last fluid plane x= 99 mm) : %8.2f Pa gauge' % p_face)
    print('  p (first bed plane x=101 mm) : %8.2f Pa gauge' % p_bed0)
    print('  p (last bed plane  x=199 mm) : %8.2f Pa gauge' % p_out)
    print('  interface jump (few Pa due to the staggered grid) : %6.2f Pa'
          % (p_bed0 - p_face))
    print('  bed drop (analytic gradient over L=0.1 m) : %8.2f Pa'
          % dp_bed)
    print('  analytic (Darcy + Forchheimer, see unMesh.control) : %.1f Pa'
          % DP_ANALYTIC)
    # |dP/dx| must match the analytic gradient (the sign is negative: the
    # pressure falls along +x, as it must for flow driven from x=0)
    ok = abs(abs(g_bed) / GRAD_ANALYTIC - 1.0) < 0.01
    print('  RESULT: %s' % ('PASS (within 1% of analytic)'
                             if ok else 'FAIL (bed gradient off by >1%)'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
