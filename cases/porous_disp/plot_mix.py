#!/usr/bin/env python3
# Plot the 2D porous mixing layer: computed T(y) column profiles vs the
# error-function similarity solution, for the dispersive and molecular-only
# runs.
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import erf

LX, LY, NX, NY = 0.04, 0.06, 80, 80
RHO, CP, K, ALP_T, U = 1.177, 1005.0, 0.026, 1.0e-3, 1.0


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
    ys = (np.arange(NY) + 0.5) * LY / NY
    prof = lambda col: np.array([T[(ix == col) & (iy == j)].mean()
                                 for j in range(NY)])
    return ys, prof


kappa = (K + RHO * CP * ALP_T * U) / (RHO * CP)      # with dispersion
kappa0 = K / (RHO * CP)                              # molecular only
ys, prof_d = load('mix_disp.vtu')
_, prof_0 = load('mix_nodisp.vtu')

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
cols = [(16, 'C0'), (40, 'C1'), (76, 'C2')]
for col, c in cols:
    x = (col + 0.5) * LX / NX
    Ta = 305 - 5 * np.array([erf((y - LY / 2) * np.sqrt(U / (4 * kappa * x)))
                             for y in ys])
    axes[0].plot(1e3 * (ys - LY / 2), prof_d(col), 'o', ms=3, color=c,
                 label='x=%.1f mm' % (1e3 * x))
    axes[0].plot(1e3 * (ys - LY / 2), Ta, '-', color=c, lw=1.1,
                 label='erf, x=%.1f mm' % (1e3 * x))
axes[0].set_title(r'disp_t = 1e-3 m: SIMPLE vs erf similarity')
axes[0].set_xlabel('y - H/2 (mm)')
axes[0].set_ylabel('T (K)')
axes[0].legend(fontsize=7)
axes[0].grid(alpha=0.3)

for col, c in cols:
    x = (col + 0.5) * LX / NX
    Ta = 305 - 5 * np.array([erf((y - LY / 2) * np.sqrt(U / (4 * kappa0 * x)))
                             for y in ys])
    axes[1].plot(1e3 * (ys - LY / 2), prof_0(col), 'o', ms=3, color=c,
                 label='x=%.1f mm' % (1e3 * x))
    axes[1].plot(1e3 * (ys - LY / 2), Ta, '-', color=c, lw=1.1,
                 label='erf, x=%.1f mm' % (1e3 * x))
axes[1].set_title(r'disp_t = 0: molecular only (under-resolved, contrast)')
axes[1].set_xlabel('y - H/2 (mm)')
axes[1].legend(fontsize=7)
axes[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig('images/mix_layer.png', dpi=140)
print('wrote images/mix_layer.png')
