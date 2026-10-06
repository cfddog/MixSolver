#!/usr/bin/env python3
"""Per-x-plane pressure/velocity dump for an uns single-solver VTU.

The uns writer (mod_uns_output.f90 vtk_write) appends one CENTROID point per
cell after the mesh nodes, and either splits every hex into 12 tetrahedra (12
sub-cells per hex, all sharing the hex centroid) or writes 5-node pyramids.
In both layouts the 5th node of each cell block is that centroid point, so the
exact per-cell centre and value are recovered without any averaging (a nodal
or sub-cell average would smear the fluid/porous interface).

This is the diagnostic behind README.md section 5 (the ~-250 Pa absolute
pressure-level offset of the single-solver reference): it prints the pmax/pmin
cell locations, the per-plane mean p/u and the fitted dP/dx over a window.

Usage:  python3 plane_profile.py unMesh.vtu [x0 [x1]]        # x in metres
"""
import re
import sys

import numpy as np


def load(path):
    txt = open(path).read()

    def raw(name):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        if m is None:
            raise SystemExit('plane_profile: no "%s" in %s' % (name, path))
        return np.fromstring(m.group(1), sep=' ')

    ncell = int(re.search(r'NumberOfCells="(\d+)"', txt).group(1))
    pts = raw('Coordinates').reshape(-1, 3)
    assert ncell > 0, 'plane_profile: empty mesh in %s' % path
    con = raw('connectivity').astype(np.int64)
    pres = raw('pressure')
    vel = raw('velocity').reshape(-1, 3)
    cen = con[4::5]                            # 5th node of every cell block
    keep = np.flatnonzero(np.concatenate(([True], cen[1:] != cen[:-1])))
    cen = cen[keep]                            # 12 sub-tets -> one hex entry
    x = pts[cen, 0]                            # centroid point ids = global
    return x, pres[keep], vel[keep, 0]


def main():
    path = sys.argv[1]
    x0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.110
    x1 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.190
    x, p, u = load(path)
    print('%s: %d cells, x = %.4f .. %.4f m'
          % (path, p.size, x.min(), x.max()))
    for tag, i in (('pmax', p.argmax()), ('pmin', p.argmin())):
        print('  %s %10.4f Pa at x=%.4f m   (u=%9.4f)'
              % (tag, p[i], x[i], u[i]))
    xr = np.round(x, 9)
    planes = np.array(sorted(set(xr)))
    mean_p = np.array([p[xr == v].mean() for v in planes])
    mean_u = np.array([u[xr == v].mean() for v in planes])
    m = (planes > x0) & (planes < x1)
    g = np.polyfit(planes[m], mean_p[m], 1)[0]
    print('  %d x-planes (dx = %.4f m); dP/dx over %.3f-%.3f m = %.1f Pa/m'
          % (planes.size, planes[1] - planes[0], x0, x1, g))
    for k in list(range(3)) + list(range(planes.size - 3, planes.size)):
        sel = xr == planes[k]
        print('    x=%8.4f mm  n=%3d  p=%10.4f  u=%9.4f'
              % (planes[k] * 1e3, int(sel.sum()), mean_p[k], mean_u[k]))


if __name__ == '__main__':
    main()
