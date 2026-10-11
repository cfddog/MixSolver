#!/usr/bin/env python3
"""Unified single-domain DBF study: compare the partitioned (legacy, explicit
fluid/porous interface treatment) results against the unified solution
(porous_model = unified: one DBF equation per cell, no explicit interface
handling, harmonic-mean face diffusion) and the Betchen/Vafai-Kim reference.

usage: compare_unified.py
Writes images/bj_unified_vs_partitioned.png and
       images/plug_unified_vs_partitioned.png, and prints a metrics table.
"""
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H = 0.01
RHO, MU = 1.177, 1.846e-5
U0 = 1.5684e-3
LX = 8 * H
NX, NY = 100, 40
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


# ---------------- BJ: fully-developed u(y) profile vs reference ------------
def bj_profile(cc, ucn):
    dx = LX / NX
    j = np.abs(cc[:, 0] - 0.95 * LX) < dx / 2
    ys = cc[j, 1] / H
    us = ucn[j]
    o = np.argsort(ys)
    return ys[o], us[o]


def bj_ref(name):
    rows = []
    with open(name) as f:
        next(f)
        for ln in f:
            if not ln.strip():
                continue
            u, y = ln.split(',')
            rows.append((float(y), float(u)))
    rows.sort(key=lambda t: t[0])
    return np.array(rows)[:, 0], np.array(rows)[:, 1]


def bj_metrics(vtu, csv):
    cc, ucn, pp_ = load(vtu)
    sy, su = bj_profile(cc, ucn)
    ry, ru = bj_ref(os.path.join(REF, csv))
    U0f = np.mean(su[sy >= 1.0])          # pure-fluid-portion mean
    yq = np.linspace(0, 2, 401)
    ui = np.interp(yq, sy, su / U0f)
    ri = np.interp(yq, ry, ru)
    l2 = np.sqrt(np.mean((ui - ri) ** 2)) / np.sqrt(np.mean(ri ** 2))
    return dict(y=yq, u=ui, l2=100 * l2, peak=ui.max(), U0f=U0f)


# ---------------- PLUG: centerline u, p + interface jitter ----------------
def plug_centerline(cc, ucn, pp):
    j = np.abs(cc[:, 1] - H / 2) < (H / 21) * 0.5
    x = cc[j, 0] / H
    uu = ucn[j] / U0
    pq = pp[j] / (RHO * U0 ** 2)
    o = np.argsort(x)
    return x[o], uu[o], pq[o]


def plug_ref(name):
    rows = []
    with open(name) as f:
        next(f)
        for ln in f:
            if not ln.strip():
                continue
            x0, v = ln.split(',')
            rows.append((float(x0), float(v)))
    rows.sort(key=lambda t: t[0])
    return np.array(rows)[:, 0], np.array(rows)[:, 1]


def plug_metrics(vtu, ucsv, pcsv):
    cc, ucn, pp = load(vtu)
    x, uu, pq = plug_centerline(cc, ucn, pp)
    rx_u, ru = plug_ref(os.path.join(REF, ucsv))
    rx_p, rp = plug_ref(os.path.join(REF, pcsv))
    # interface-neighbourhood odd-even jitter
    win = (x > 2.3) & (x < 5.8)
    ii = np.where(win)[0][1:-1]
    alt = np.abs(uu[ii] - 0.5 * (uu[ii - 1] + uu[ii + 1])).max()
    ui = np.interp(rx_u, x, uu)
    pi = np.interp(rx_p, x, pq)
    ku, kp = rx_u >= x[0], rx_p >= x[0]
    m = rx_u[ku] >= 2.0
    l2u = 100 * np.sqrt(np.mean((ui[ku][m] - ru[ku][m]) ** 2)) / np.abs(ru[ku][m]).max()
    l2p = 100 * np.sqrt(np.mean((pi[kp] - rp[kp]) ** 2)) / (np.abs(rp).max() + 1e-30)
    return dict(x=x, u=uu, p=pq, alt=alt, l2u=l2u, l2p=l2p)


# ============================ report ============================
print('=' * 78)
print('BJ (parallel interface, fully-developed u(y), L2 on shape)')
print('-' * 78)
print('%-10s %-16s %-16s' % ('case', 'partitioned', 'unified'))
fig_bj, axb = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for k, (tag, csv) in enumerate([('bj_dae2', 'BJ_1.csv'), ('bj_dae3', 'BJ_2.csv')]):
    mp = bj_metrics('%s.vtu' % tag, csv)
    mu = bj_metrics('%s_uni.vtu' % tag, csv)
    print('%-10s L2=%6.3f%%      L2=%6.3f%%      (U0f/U0in p=%.3f u=%.3f)'
          % (tag, mp['l2'], mu['l2'], mp['U0f'] / U0, mu['U0f'] / U0))
    ry, ru = bj_ref(os.path.join(REF, csv))
    ax = axb[k]
    ax.plot(mp['u'], mp['y'], '-o', ms=3, lw=1.4, label='partitioned')
    ax.plot(mu['u'], mu['y'], '-s', ms=3, lw=1.4, label='unified')
    ax.plot(ru, ry, '--', lw=1.8, color='C1', label='Betchen ref')
    ax.set_title(tag)
    ax.set_xlabel('u / U0')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
axb[0].set_ylabel('y / H')
fig_bj.tight_layout()
fig_bj.savefig('images/bj_unified_vs_partitioned.png', dpi=140)
print('wrote images/bj_unified_vs_partitioned.png')

print()
print('=' * 78)
print('PLUG (perpendicular interface: jitter alt_amp, centerline L2)')
print('-' * 78)
print('%-10s %-22s %-22s' % ('case', 'partitioned', 'unified'))
print('%-10s %-22s %-22s' % ('', 'alt / uL2% / pL2%', 'alt / uL2% / pL2%'))
fig_pl, axp = plt.subplots(1, 2, figsize=(12, 4.6))
for k, (tag, ucsv, pcsv, col) in enumerate(
        [('plug_dae2', 'plug_1_u.csv', 'plug_1_p.csv', 'C0'),
         ('plug_dae3', 'plug_2_u.csv', 'plug_2_p.csv', 'C3')]):
    mp = plug_metrics('%s.vtu' % tag, ucsv, pcsv)
    mu = plug_metrics('%s_uni.vtu' % tag, ucsv, pcsv)
    print('%-10s %.4f / %6.3f / %6.3f   %.4f / %6.3f / %6.3f'
          % (tag, mp['alt'], mp['l2u'], mp['l2p'],
             mu['alt'], mu['l2u'], mu['l2p']))
    rx_u, ru = plug_ref(os.path.join(REF, ucsv))
    rx_p, rp = plug_ref(os.path.join(REF, pcsv))
    axp[0].plot(mp['x'], mp['u'], '-o', ms=4, lw=1.6, color=col,
                fillstyle='none', label='partitioned ' + tag)
    axp[0].plot(mu['x'], mu['u'], '-s', ms=4, lw=1.6, color=col,
                fillstyle='none', label='unified ' + tag)
    axp[0].plot(rx_u, ru, '--', lw=1.8, color=col, label='ref ' + tag)
    axp[1].plot(mp['x'], mp['p'], '-o', ms=4, lw=1.6, color=col,
                fillstyle='none', label='partitioned ' + tag)
    axp[1].plot(mu['x'], mu['p'], '-s', ms=4, lw=1.6, color=col,
                fillstyle='none', label='unified ' + tag)
    axp[1].plot(rx_p, rp, '--', lw=1.8, color=col, label='ref ' + tag)
axp[0].set_xlabel('x / H')
axp[0].set_ylabel('u / U (centerline)')
axp[0].set_title('Centerline velocity')
axp[0].grid(alpha=0.3)
axp[0].legend(fontsize=7, ncol=2)
axp[1].set_xlabel('x / H')
axp[1].set_ylabel('p / (rho U^2)')
axp[1].set_title('Centerline pressure')
axp[1].grid(alpha=0.3)
axp[1].legend(fontsize=7, ncol=2)
fig_pl.tight_layout()
fig_pl.savefig('images/plug_unified_vs_partitioned.png', dpi=140)
print('wrote images/plug_unified_vs_partitioned.png')
