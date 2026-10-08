#!/usr/bin/env python3
# Betchen PLUG post-processing: velocity and cross-section-mean pressure along
# x.  Tests the paper's claim of a LINEAR pressure profile in each of the
# three segments (fluid/porous/fluid) with a pressure-gradient discontinuity
# at the two interfaces.  Reports per-segment linear-fit R^2 and slopes.
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H = 0.01


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
    N = len(off)                     # 12 tets per hex
    c = np.zeros((N, 3))
    for i in range(N):
        ids = conn[st[i]:off[i]]
        if ids[0] == 4 and len(ids) == 5:
            ids = ids[1:]
        c[i] = pts[ids].mean(0)
    nhex = N // 12
    cc = c.reshape(nhex, 12, 3).mean(1)
    uU = U[:N].reshape(nhex, 12, 3).mean(1)[:, 0]
    pp = p[:N].reshape(nhex, 12).mean(1)
    return cc, uU, pp


def column_means(cc, uU, pp, nx):
    x0 = cc[:, 0].min(); x1 = cc[:, 0].max()
    w = (x1 - x0) / nx
    ix = np.floor((cc[:, 0] - x0) / w).astype(int)
    ix = np.clip(ix, 0, nx - 1)
    xm = np.array([cc[ix == i, 0].mean() if (ix == i).any() else np.nan
                   for i in range(nx)])
    um = np.array([uU[ix == i].mean() if (ix == i).any() else np.nan
                   for i in range(nx)])
    pm = np.array([pp[ix == i].mean() if (ix == i).any() else np.nan
                   for i in range(nx)])
    return xm, um, pm


def fit(x, y, xa, xb):
    m = (x >= xa) & (x <= xb) & np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan, 0.0, m.sum()
    A = np.vstack([x[m], np.ones(len(x[m]))]).T
    (s, c0), _, _, _ = np.linalg.lstsq(A, y[m], rcond=None)
    r2 = 1.0 - np.sum((y[m] - (s * x[m] + c0)) ** 2) / \
         (np.sum((y[m] - y[m].mean()) ** 2) + 1e-30)
    return s, r2, m.sum()


cases = [
    ('plug_dae2.vtu', 'PLUG Da=1e-2 Re=1', (3, 2, 3)),
    ('plug_dae3.vtu', 'PLUG Da=1e-3 Re=1', (3, 2, 3)),
    ('plug_hir.vtu', 'PLUG Da=1e-2 Re=1000', (5, 5, 50)),
]

fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
for fn, lab, seg in cases:
    xhi = sum(s * H for s in seg)
    cc, uU, pp = load(fn)
    xm, um, pm = column_means(cc, uU, pp, 150)
    axs[0].plot(1e2 * xm, pm, '-o', ms=2, lw=1.2, label=lab)
    axs[1].plot(1e2 * xm, um, '-o', ms=2, lw=1.2, label=lab)
    # per-segment linear fit of the pressure profile
    x1 = seg[0] * H; x2 = x1 + seg[1] * H
    print('== %s  porous x in [%.3f, %.3f] m' % (lab, x1, x2))
    for si, (xa, xb, nm) in enumerate(
            [(0.0, x1, 'seg1 fluid'), (x1, x2, 'seg2 porous'),
             (x2, xhi, 'seg3 fluid')]):
        s, r2, n = fit(xm, pm, xa, xb)
        print('   %-11s x[%.3f,%.3f]  n=%-3d  slope=%10.4e Pa/m  R^2=%.4f'
              % (nm, xa, xb, n, s, r2))
    # mark porous segment
    axs[0].axvspan(x1 * 1e2, x2 * 1e2, color='gray', alpha=0.15)
    axs[1].axvspan(x1 * 1e2, x2 * 1e2, color='gray', alpha=0.15)

axs[0].set_xlabel('x (cm)'); axs[0].set_ylabel('p mean (Pa)')
axs[0].set_title('Pressure along x')
axs[0].grid(alpha=0.3); axs[0].legend(fontsize=7)
axs[1].set_xlabel('x (cm)'); axs[1].set_ylabel('u mean (m/s)')
axs[1].set_title('Velocity along x')
axs[1].grid(alpha=0.3); axs[1].legend(fontsize=7)
fig.tight_layout()
fig.savefig('images/plug_pressure.png', dpi=140)
print('wrote images/plug_pressure.png')