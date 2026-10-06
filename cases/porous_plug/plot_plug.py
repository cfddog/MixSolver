#!/usr/bin/env python3
# Plot the 1D porous-plug results: pressure column means vs the analytical
# Darcy / Darcy-Forchheimer solutions, and a normalised velocity profile.
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LX, NX = 0.1, 40
RHO, MU, EPS, K = 1.177, 1.846e-5, 0.4, 1.0e-8


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
    U = arr('velocity', cd_).reshape(-1, 3)
    p = arr('pressure', cd_)
    pts = arr('Coordinates', re.search(r'<Points>(.*?)</Points>', s, re.S)
              .group(1)).reshape(-1, 3)
    st = np.concatenate(([0], off[:-1]))
    xc = np.zeros(len(off))
    for i in range(len(off)):                 # tets: strip leading n-verts token
        ids = conn[st[i]:off[i]]
        if ids[0] == 4 and len(ids) == 5:
            ids = ids[1:]
        xc[i] = pts[ids].mean(0)[0]
    ix = np.round((xc - LX / NX / 2) / (LX / NX)).astype(int)
    xm = np.array([(i + 0.5) * LX / NX for i in range(NX)])
    um = np.array([U[ix == i, 0].mean() for i in range(NX)])
    pm = np.array([p[ix == i].mean() for i in range(NX)])
    return xm, um, pm


cases = [
    ('plug_darcy.vtu', 0.1, 0.0, 'C0', 'Darcy only ($u_{in}=0.1$ m/s)'),
    ('plug_forch.vtu', 1.0, 500.0, 'C1', 'Darcy-Forchheimer ($u_{in}=1$ m/s)'),
]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
for fn, u_in, cf, c, lab in cases:
    x, u, p = load(fn)
    ax1.plot(1e3 * x, p, 'o', ms=4, color=c, label=lab + ' (SIMPLE)')
    dpdx = -(MU / EPS / K) * u_in - RHO * cf * u_in**2
    p_th = dpdx * (x - LX)         # zero pressure at the outlet face x=L
    ax1.plot(1e3 * x, p_th, '-', color=c, lw=1.2,
             label=lab.split('(')[0] + 'theory')
    ax2.plot(1e3 * x, u / u_in, 'o', ms=4, color=c, label=lab)

ax1.set_xlabel('x (mm)')
ax1.set_ylabel('p gauge (Pa)')
ax1.set_title('Pressure: SIMPLE column means vs analytical')
ax1.legend(fontsize=8)
ax1.grid(alpha=0.3)
ax2.axhline(1.0, color='k', ls='--', lw=1)
ax2.set_xlabel('x (mm)')
ax2.set_ylabel(r'$u/u_{in}$')
ax2.set_ylim(0.85, 1.15)
ax2.set_title('Velocity (inlet/outlet column artifacts excluded interior)')
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('images/plug_pressure.png', dpi=140)
print('wrote images/plug_pressure.png')
