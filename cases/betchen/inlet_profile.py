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
# y/H = 1/(2*NY) and 1-1/(2*NY); int(u dy)/(U0 H) over that range is printed
# as a rough check that the profile mean is U0 (it approaches 1 for a
# developed/parabolic profile, and is < 1 right at a flat inlet because the
# wall cells are already being dragged).
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
        print('   x/H=%.3f  y/H : %s'
              % (xc, ' '.join('%5.2f' % v for v in xi)))
        print('             ana : %s'
              % ' '.join('%5.3f' % v for v in 6 * xi * (1 - xi)))
        print('             our : %s'
              % ' '.join('%5.3f' % v for v in uu[o]))
        print('             int(u dy)/(U0 H) = %.5f   (target 1.00000)'
              % (np.trapezoid(uu[o], yy[o])))
