#!/usr/bin/env python3
# compact comparison of Betchen PLUG runs: jitter amplitude + pressure level
import re
import numpy as np

H = 0.01
RHO, MU = 1.177, 1.846e-5
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
    p = arr('pressure', cd_)
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
    ucn = U[:N].reshape(nhex, 12, 3).mean(1)[:, 0]
    pp = p[:N].reshape(nhex, 12).mean(1)
    return cc, ucn, pp


def report(fn, lab, u0=None):
    u0 = U0 if u0 is None else u0
    cc, u, p = load(fn)
    xs = cc[:, 0] / H
    j = np.abs(cc[:, 1] - H / 2) < (H / 21) * 0.5
    x, uu, pq = xs[j], u[j] / u0, p[j] / (RHO * u0 ** 2)
    o = np.argsort(x)
    x, uu, pq = x[o], uu[o], pq[o]
    # odd-even (cell-to-cell) amplitude in the interface neighbourhood
    win = (x > 2.3) & (x < 5.8)
    ii = np.where(win)[0][1:-1]
    alt = np.abs(uu[ii] - 0.5 * (uu[ii - 1] + uu[ii + 1]))
    # plug gradient (interior of the 2H plug at x/H in [3,5])
    m = (x > 3.25) & (x < 4.75)
    g = np.polyfit(x[m], pq[m], 1)[0]
    k = np.argsort(np.abs(x - 4.05))[0]
    print('%-8s  alt_amp=%.4f  (u range in win %.3f..%.3f)  '
          'plug dpdx=%8.1f /H  p_inlet=%8.1f  p(x=4.05)=%8.1f'
          % (lab, alt.max(), uu[win].min(), uu[win].max(), -g, pq[0], pq[k]))
    print('         centerline u/U x/H 2.3..5.8: '
          + ' '.join('%.3f' % v for v in uu[win]))
    print('         centerline p     x/H 2.3..5.8: '
          + ' '.join('%.1f' % v for v in pq[win]))


import sys
for a in sys.argv[1:]:
    parts = a.split('=')
    fn, lab = parts[0], parts[1]
    u0 = float(parts[2]) if len(parts) > 2 else None
    report(fn, lab, u0)
