#!/usr/bin/env python3
# Compare the Betchen BJ trial-run fully-developed profile against the provided
# reference (BJ_1.csv = Da=1e-2, BJ_2.csv = Da=1e-3).  Both are u(y) curves at
# the fully-developed section; the reference uses a normalizing velocity ~2x
# U0 (convention), so we fit a constant scale to align peaks, then report the
# residual on the SHAPE.  Writes images/bj_compare_ref.png.
import re, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H = 0.01
LX = 8 * H
NX, NY = 100, 40
U0 = 1.5684e-3
REF = '/mnt/c/temp/validate_case'


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
    uU = U[:N].reshape(nhex, 12, 3).mean(1)[:, 0]
    return cc, uU


def profile(cc, uU):
    dx = LX / NX
    x0 = 0.95 * LX                       # fully-developed (peak/area stable)
    j = np.abs(cc[:, 0] - x0) < dx / 2
    ys = cc[j, 1] / H
    us = uU[j]
    o = np.argsort(ys)
    return ys[o], us[o]


def get_ref(name):
    rows = []
    with open(name) as f:
        next(f)
        for ln in f:
            if not ln.strip():
                continue
            u, y = ln.split(',')
            rows.append((float(y), float(u)))
    rows.sort(key=lambda t: t[0])            # sort by y
    return np.array(rows)[:, 0], np.array(rows)[:, 1]   # (y, u)
cases = [('bj_dae2.vtu', os.path.join(REF, 'BJ_1.csv'), 'Da=1e-2'),
         ('bj_dae3.vtu', os.path.join(REF, 'BJ_2.csv'), 'Da=1e-3')]

fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, (vtu, csv, lab) in zip(axs, cases):
    cc, uU = load(vtu)
    sy, su = profile(cc, uU)             # su in m/s, y in units of H
    ry, ru = get_ref(csv)                # reference u/U0 (U0 = fluid-portion mean)
    # U0 = fully-developed average velocity in the PURE FLUID portion (y>H)
    fl = sy >= 1.0
    U0f = np.mean(su[fl])
    yq = np.linspace(0, 2, 401)
    ui = np.interp(yq, sy, su / U0f)     # mine, normalized by the same U0
    ri = np.interp(yq, ry, ru)
    resid = ui - ri
    rms = np.sqrt(np.mean(ri ** 2))
    l2 = np.sqrt(np.mean(resid ** 2))
    print('== %s  U0_fluid=%.4e m/s (%.3f * U0_inlet): L2 residual=%.4f '
          '(%.2f%% of ref RMS %.4f)' % (lab, U0f, U0f / U0, l2, 100 * l2 / rms,
                                        rms))
    print('   my-fluid-mean-peak=%.3f vs ref-peak=%.3f' % (ui.max(), ru.max()))
    ax.plot(ui, yq, '-o', ms=3, lw=1.4, label='SIMPLE (trial)')
    ax.plot(ru, ry, '--s', ms=4, color='C1', label='Betchen ref ' + lab)
    ax.set_title(lab)
    ax.set_xlabel('u / U0')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
axs[0].set_ylabel('y / H')
fig.tight_layout()
fig.savefig('images/bj_compare_ref.png', dpi=140)
print('wrote images/bj_compare_ref.png')