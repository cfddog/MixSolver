#!/usr/bin/env python3
"""Geometry-free full-field diff of two VTUs (BJ meshes have a different
aspect ratio than the PLUG ones, so cmp.py/diff2.py's H-based centerline
selection does not apply).

usage: fdiff.py a.vtu b.vtu [label]
"""
import re
import sys

import numpy as np


def load(fn):
    s = open(fn).read()

    def arr(name, text):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii">(.*?)'
                      r'</DataArray>' % re.escape(name), text, re.S)
        return np.fromstring(m.group(1), sep=' ')

    cd = re.search(r'<CellData>(.*?)</CellData>', s, re.S).group(1)
    u = arr('velocity', cd).reshape(-1, 3)
    p = arr('pressure', cd)
    return u, p


a, b = sys.argv[1], sys.argv[2]
lab = sys.argv[3] if len(sys.argv) > 3 else ''
ua, pa = load(a)
ub, pb = load(b)
n = min(len(ua), len(ub))
du = np.abs(ua[:n] - ub[:n])
dp = np.abs(pa[:n] - pb[:n])
umax = np.abs(ua[:n]).max()
pmax = np.abs(pa[:n]).max()
print('%-22s  cells=%d' % (lab or (a + ' vs ' + b), n))
print('   full field: max|du| = %.4e (%.4f%% of u_max=%.4e)   '
      'L2 = %.4e  max|dp| = %.4e (%.4f%% of p_max=%.4e)  L2 = %.4e'
      % (du.max(), 100 * du.max() / umax, umax,
         np.sqrt((du ** 2).mean()),
         dp.max(), 100 * dp.max() / pmax, pmax,
         np.sqrt((dp ** 2).mean())))
