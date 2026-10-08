#!/usr/bin/env python3
"""Pressure-gradient fit for the 1-D porous-plug cases (cases/porous_plug).

Reads a VTU (ASCII), reconstructs cell centres from the node coordinates,
averages each x-column, fits p(x) over the interior columns (col 3..36, i.e.
x/L = 0.075..0.90) and prints the fitted dp/dx together with the analytic
value for the *Nield-Bejan / Fluent* convention  dp/dx = -(mu/K) u  plus the
Forchheimer term -(rho*inertial)|u|u.
"""
import re
import sys

import numpy as np

RHO = 1.177
MU = 1.846e-5
L = 0.1


def load(fn):
    s = open(fn).read()

    def arr(name, text):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii">(.*?)'
                      r'</DataArray>' % re.escape(name), text, re.S)
        return np.fromstring(m.group(1), sep=' ')

    cell = re.search(r'<Cells>(.*?)</Cells>', s, re.S).group(1)
    conn = arr('connectivity', cell).astype(int)
    off = arr('offsets', cell).astype(int)
    cd = re.search(r'<CellData>(.*?)</CellData>', s, re.S).group(1)
    p = arr('pressure', cd)
    U = arr('velocity', cd).reshape(-1, 3)
    pts = arr('Coordinates', re.search(r'<Points>(.*?)</Points>', s, re.S)
              .group(1)).reshape(-1, 3)
    st = np.concatenate(([0], off[:-1]))
    n = len(off)
    xc = np.zeros(n)
    for i in range(n):
        ids = conn[st[i]:off[i]]
        if ids[0] == 4 and len(ids) == 5:
            ids = ids[1:]
        xc[i] = pts[ids, 0].mean()
    return xc, p, U[:n]


def fit(fn, u_in, K, eps, inertial=0.0):
    xc, p, U = load(fn)
    # column index: the mesh is 40 cells in x, uniform
    icol = np.round(xc / L * 40.0).astype(int)
    keep = (icol >= 3) & (icol <= 36)
    x = np.array([xc[keep & (icol == c)].mean() for c in np.unique(icol[keep])])
    pv = np.array([p[keep & (icol == c)].mean() for c in np.unique(icol[keep])])
    uu = np.array([U[keep & (icol == c), 0].mean()
                   for c in np.unique(icol[keep])])
    a, b = np.polyfit(x, pv, 1)
    a_ana = -(MU / K) * u_in - RHO * inertial * abs(u_in) * u_in
    print('%-16s u_in=%.4g  fitted dp/dx = %10.3f   analytic(%s) = %10.3f  '
          'err=%6.3f%%   u max dev = %.4f%%   dp = %8.4f Pa'
          % (fn.split('/')[-1], u_in, a,
             'mu/K + Forch' if inertial else 'mu/K', a_ana,
             100.0 * (a - a_ana) / a_ana,
             100.0 * np.max(np.abs(uu - u_in)) / u_in,
             -a * L))


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        name = fn.split('/')[-1]
        if 'forch' in name:
            fit(fn, 1.0, 1.0e-8, 0.4, inertial=500.0)
        elif 'aniso2' in name:
            fit(fn, 0.1, 5.0e-9, 0.4)
        else:
            fit(fn, 0.1, 1.0e-8, 0.4)
