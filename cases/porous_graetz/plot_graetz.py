#!/usr/bin/env python3
# Graetz validation plots:
#   left/mid : theta(eta) column profiles vs the plug-flow Graetz series
#   right    : bulk theta vs x* for both runs and the series (collapse check)
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import pi

LX, LY, NX, NY = 0.45, 0.005, 180, 40
A = LY
RHO, CP, K = 1.177, 1005.0, 0.026


def load(fn):
    s = open(fn).read()

    def arr(name, text):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii">(.*?)</DataArray>'
                      % re.escape(name), text, re.S)
        return np.fromstring(m.group(1), sep=' ')

    sec = re.search(r'<Cells>(.*?)</Cells>', s, re.S).group(1)
    conn = arr('connectivity', sec).astype(int)
    off = arr('offsets', sec).astype(int)
    cd_ = re.search(r'<CellData>(.*?)</CellData>', s, re.S).group(1)
    T = arr('temperature', cd_)
    pts = arr('Coordinates', re.search(r'<Points>(.*?)</Points>', s, re.S)
              .group(1)).reshape(-1, 3)
    st = np.concatenate(([0], off[:-1]))
    xyc = np.zeros((len(off), 2))
    for i in range(len(off)):                 # tets: strip leading n-verts token
        ids = conn[st[i]:off[i]]
        if ids[0] == 4 and len(ids) == 5:
            ids = ids[1:]
        xyc[i] = pts[ids].mean(0)[:2]
    ix = np.round((xyc[:, 0] - LX / NX / 2) / (LX / NX)).astype(int)
    iy = np.round((xyc[:, 1] - LY / NY / 2) / (LY / NY)).astype(int)
    return ix, iy, T


lam = np.array([(n - 0.5) * pi for n in range(1, 80)])
sgn = np.array([(-1) ** (n + 1) for n in range(1, 80)])
Cn = 2 * sgn / lam
eta_f = np.linspace(0, 1, 200)


def tprof(xs, eta):
    return np.sum(Cn * np.cos(lam * eta) * np.exp(-lam ** 2 * xs))


def tbulk(xs):
    return np.sum(2.0 / lam ** 2 * np.exp(-lam ** 2 * xs))


eta = (np.arange(NY) + 0.5) / NY
fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
runs = [('graetz_mol.vtu', 0.1, 0.0, axes[0], 'molecular, Pe=22.7'),
        ('graetz_disp.vtu', 1.0, 2.5e-4, axes[1], 'disp_t=2.5e-4 m, Pe=18.4')]
for fn, U, alpt, ax, title in runs:
    ix, iy, T = load(fn)
    kap = (K + RHO * CP * alpt * U) / (RHO * CP)
    for xs, c in [(0.2, 'C0'), (0.5, 'C1'), (1.0, 'C2'), (2.0, 'C3')]:
        col = int(round(xs * U * A ** 2 / (kap * LX / NX) - 0.5))
        xs_a = kap * ((col + 0.5) * LX / NX) / (U * A ** 2)
        Tc = np.array([T[(ix == col) & (iy == j)].mean() for j in range(NY)])
        th = (Tc - 310.0) / -10.0
        ax.plot(th, eta, 'o', ms=3, color=c)
        ax.plot([tprof(xs_a, e) for e in eta_f], eta_f, '-', color=c, lw=1,
                label='x*=%.2f' % xs_a)
    ax.set_title(title)
    ax.set_xlabel(r'$\theta=(T-T_w)/(T_0-T_w)$')
    ax.set_ylabel(r'$\eta=y/a$')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

ax = axes[2]
xx = np.linspace(0.02, 4.2, 200)
ax.plot(xx, [tbulk(v) for v in xx], 'k-', lw=1.2, label='Graetz series')
for fn, U, alpt, m in [('graetz_mol.vtu', 0.1, 0.0, 's'),
                        ('graetz_disp.vtu', 1.0, 2.5e-4, '^')]:
    ix, iy, T = load(fn)
    kap = (K + RHO * CP * alpt * U) / (RHO * CP)
    cols = np.arange(3, NX - 2, 4)
    xs = kap * (cols + 0.5) * LX / NX / (U * A ** 2)
    thb = np.array([T[ix == c].mean() for c in cols])  # uniform u -> column mean
    thb = (thb - 310.0) / -10.0
    ax.plot(xs, thb, m, ms=4,
            label=('molecular' if alpt == 0 else 'with dispersion'))
ax.set_xlabel(r'$x^*=\kappa_{eff}x/(Ua^2)$')
ax.set_ylabel(r'bulk $\theta_b$')
ax.set_title('Bulk temperature collapse')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('images/graetz.png', dpi=140)
print('wrote images/graetz.png')
