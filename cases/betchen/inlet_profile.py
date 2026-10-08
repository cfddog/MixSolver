#!/usr/bin/env python3
# Print the first few x-columns of the z-center plane as u(y)/U0 profiles, to
# check the inlet boundary condition shape directly against the analytic
# fully-developed parabola 6*xi*(1-xi) (xi = y/H).  Use it to compare a
# 'velocity-inlet' (flat) run with a 'velocity-inlet-parabolic' one.
#
# usage: python3 inlet_profile.py <vtu>[=<label>][=<U0>] [more...]
#   e.g. python3 inlet_profile.py /tmp/par/plug_dae3.vtu=PAR
#        python3 inlet_profile.py plug_dae3.vtu=PAR 1.5684e-3
#
# Note: velocities are cell-centred, so the outermost samples are at
# y/H = 1/(2*NY) and 1-1/(2*NY).  int(u dy)/(U0 H) is printed together with the
# SAME discrete integral of the analytic parabola sampled at the same y's
# ('disc' below): the continuous target 1.0 is NOT reachable from cell-centred
# data because the wall margins y/H<1/(2NY) and >1-1/(2NY) are missing (they
# cost ~0.55% at NY=21, i.e. disc=0.99449).  Compare 'our' with 'disc', not 1.0.
#
# The FIRST column is the inlet cell: the profile is imposed exactly on the face
# at x=0, but the first cell centre sits dx/2 downstream and its value carries an
# O(dx) inlet-cell error (centre u/U0 = 1.483 at dx=0.15H, 1.497 at dx=0.05H),
# so do not read the first column as "the BC is wrong" -- refine x instead.
import re
import sys
import numpy as np

H = 0.01
U0 = 1.5684e-3


def load(fn):
    s = open(fn).read()

    def arr(name, text):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii">(.*?)'
                      r'</DataArray>' % re.escape(name), text, re.S)
        return np.fromstring(m.group(1), sep=' ')

    sec = re.search(r'<Cells>(.*?)</Cells>', s, re.S).group(1)
    conn = arr('connectivity', sec).astype(int)
    off = arr('offsets', sec).astype(int)
    cd_ = re.search(r'<CellData>(.*?)</CellData>', s, re.S).group(1)
    U = arr('velocity', cd_).reshape(-1, 3)
    pts = arr('Coordinates', re.search(r'<Points>(.*?)</Points>', s, re.S)
              .group(1)).reshape(-1, 3)
    st = np.concatenate(([0], off[:-1]))
    N = len(off)
    c = np.zeros((N, 3))
    for i in range(N):
        ids = conn[st[i]:off[i]]
        if ids[0] == 4 and len(ids) == 5:
            ids = ids[1:]
        c[i] = pts[ids].mean(0)
    nhex = N // 12
    cc = c.reshape(nhex, 12, 3).mean(1)
    return cc, U[:N].reshape(nhex, 12, 3).mean(1)[:, 0]


for a in sys.argv[1:]:
    parts = a.split('=')
    fn, lab = parts[0], parts[1]
    u0 = float(parts[2]) if len(parts) > 2 else U0
    cc, u = load(fn)
    xs = cc[:, 0] / H
    jz = np.abs(cc[:, 2] - cc[:, 2].mean()) < 1e-9      # z-center plane
    cols = np.unique(np.round(xs[jz], 3))
    print('== %s  (first 4 x-columns of the z-center plane, u/U0)' % lab)
    for xc in cols[:4]:
        m = jz & (np.abs(xs - xc) < 1e-6)
        yy = cc[m, 1] / H
        uu = u[m] / u0
        o = np.argsort(yy)
        xi = yy[o]
        ana = 6 * xi * (1 - xi)
        print('   x/H=%.3f  y/H : %s'
              % (xc, ' '.join('%5.2f' % v for v in xi)))
        print('             ana : %s'
              % ' '.join('%5.3f' % v for v in ana))
        print('             our : %s'
              % ' '.join('%5.3f' % v for v in uu[o]))
        print('             max|our-ana| = %.4f   int(u dy)/(U0 H) = %.5f '
              '(disc %.5f)'
              % (np.abs(uu[o] - ana).max(), np.trapezoid(uu[o], yy[o]),
                 np.trapezoid(ana, xi)))
