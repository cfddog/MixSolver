#!/usr/bin/env python3
# Betchen HT post-processing: report average Tf/Ts, bottom-wall heat flux per
# unit depth q' [W/m] and wall normal temperature gradient, for both cases.
# Compare q' with the paper target 315 W/m and Nu ~ 5.1.
import re
import numpy as np

RHO, MU = 1.177, 1.846e-5
CP, KF = 1005.0, 0.026
H_ht = 0.045
EPS = 0.9118
K_SE, K_FE = 6.46, 0.0237          # solid/fluid effective conductivities


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
    Tf = arr('temperature', cd_)
    Ts = arr('temperature_solid', cd_)
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
    tf = Tf[:N].reshape(nhex, 12).mean(1)
    ts = Ts[:N].reshape(nhex, 12).mean(1)
    return cc, tf, ts


def wall_flux(cc, tf, ts, dx, dy):
    # bottom row (y ~ dy/2), per unit depth. q'_ [W/m] = sum over the bottom
    # columns of k_eff*(T_w - T_cell)/(dy/2) * dx  (wall at y=0, cell centre at
    # y=dy/2).  T_w = 310 K constant heated wall.
    bot = cc[:, 1] < dy
    qf = K_FE * (310.0 - tf[bot]) / (dy / 2) * dx
    qs = K_SE * (310.0 - ts[bot]) / (dy / 2) * dx
    return qf.sum(), qs.sum()


for fn, lab in [('ht_a.vtu', 'HT case A (block)'),
                ('ht_b.vtu', 'HT case B (gap+block)')]:
    cc, tf, ts = load(fn)
    dx = cc[:, 0].max() / 80.0
    dy = H_ht / 70.0
    print('== %s  (dx=%.4f dy=%.4f mm)' % (lab, dx * 1e3, dy * 1e3))
    print('   Tf  range [%.2f, %.2f]  mean %.2f'
          % (tf.min(), tf.max(), tf.mean()))
    print('   Ts  range [%.2f, %.2f]  mean %.2f'
          % (ts.min(), ts.max(), ts.mean()))
    print('   (Tf-Ts) hot-bottom row mean=%.3f K' %
          (tf[cc[:, 1] < dy].mean() - ts[cc[:, 1] < dy].mean()))
    qf, qs = wall_flux(cc, tf, ts, dx, dy)
    print('   bottom heat flux per depth: qf\'=%.1f W/m  qs\'=%.1f W/m  '
          'total=%.1f W/m   (paper target 315 W/m)' % (qf, qs, qf + qs))