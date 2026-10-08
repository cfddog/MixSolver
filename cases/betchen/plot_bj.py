#!/usr/bin/env python3
# Betchen BJ case post-processing: extract u(y) fully-developed profile at a
# downstream station x ~ 0.8*L, split fluid/porous at y=H, report u(y)/U0 and
# the slip velocity / Darcy velocity.  Console + PNG profile.
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H = 0.01          # half channel height
LX = 8 * H        # channel length
NX, NY = 100, 40


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
    N = len(off)                     # number of CARDS (12 tets per hex)
    c = np.zeros((N, 3))             # tet centroids
    for i in range(N):
        ids = conn[st[i]:off[i]]
        if ids[0] == 4 and len(ids) == 5:
            ids = ids[1:]
        c[i] = pts[ids].mean(0)
    nhex = N // 12
    cc = c.reshape(nhex, 12, 3).mean(1)      # hex centroids
    uU = U[:N].reshape(nhex, 12, 3).mean(1)[:, 0]
    pp = p[:N].reshape(nhex, 12).mean(1)
    return cc, uU, pp


def profile(cc, uU, pp, x0, nbin=200):
    dx = LX / NX
    j = np.abs(cc[:, 0] - x0) < dx / 2       # one cell column
    ys = cc[j, 1]
    us = uU[j]
    ps = pp[j]
    order = np.argsort(ys)
    return ys[order], us[order], ps[order]


cases = [
    ('bj_dae2.vtu', 1.0e-6, 'Da=1e-2 (K=1e-6)'),
    ('bj_dae3.vtu', 1.0e-7, 'Da=1e-3 (K=1e-7)'),
]

fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
for fn, K, lab in cases:
    cc, uU, pp = load(fn)
    x0 = 0.85 * LX
    y, u, p = profile(cc, uU, pp, x0)
    u0 = 1.5684e-3
    print('== %s  (x0=%.3f m)' % (lab, x0))
    sl = (y >= H) & (y <= 0.9 * (2 * H))
    po = y < H
    print('   fluid column mean u/U0    = %.4f  (%.2f cells)'
          % (u[sl].mean() / u0, sl.sum()))
    print('   porous column mean u/U0   = %.4f  (%.2f cells)'
          % (u[po].mean() / u0, po.sum()))
    # interface value (y just below H) and top (y~2H)
    yH = np.argmin(np.abs(y - H))
    yH1 = np.argmin(np.abs(y - (H * 1.02)))
    print('   u(interface ~ y=H) /U0    = %.4f' % (u[yH] / u0))
    axs[0].plot(1e3 * y, u / u0, '-', lw=1.8,
                label=lab + '  (u/Yfluid=%.3f, u|_bed=%.3f)'
                % (u[sl].mean() / u0, u[po].mean() / u0))
    axs[0].axhline(1.0, color='k', ls='--', lw=0.8)
    axs[1].plot(1e3 * y, p, '-', lw=1.8, label=lab)

axs[0].axhline(H * 1e3, color='k', ls=':', lw=0.8)
axs[0].set_xlabel('y (mm)')
axs[0].set_ylabel('u / U0')
axs[0].set_title('Fully-developed u(y)/U0 at ~0.85L (fluid above / porous below)')
axs[0].grid(alpha=0.3)
axs[0].legend(fontsize=8)
axs[1].set_xlabel('y (mm)')
axs[1].set_ylabel('p gauge (Pa)')
axs[1].set_title('Pressure across height at ~0.85L')
axs[1].grid(alpha=0.3)
axs[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig('images/bj_profiles.png', dpi=140)
print('wrote images/bj_profiles.png')