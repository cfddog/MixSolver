#!/usr/bin/env python3
# Betchen PLUG: our VTU centerline vs the reference CSVs (Vafai-Kim analytical).
# usage: vs_ref.py <vtu>=<ucsv>=<pcsv>=<label> [more quads...]
#   e.g.  python3 vs_ref.py plug_dae3.vtu=plug_2_u.csv=plug_2_p.csv=post
# The CSV names are resolved against $REF (default /mnt/c/temp/validate_case).
import os
import re
import sys
import numpy as np

H = 0.01
RHO, MU = 1.177, 1.846e-5
U0 = 1.5684e-3
REF = '/mnt/c/temp/validate_case'


def arr(name, text):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii">(.*?)'
                  r'</DataArray>' % re.escape(name), text, re.S)
    return np.fromstring(m.group(1), sep=' ')


def centerline(fn):
    s = open(fn).read()
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
    j = np.abs(cc[:, 1] - H / 2) < (H / 21) * 0.5
    x = cc[j, 0] / H
    o = np.argsort(x)
    return x[o], ucn[j][o] / U0, pp[j][o] / (RHO * U0 ** 2)


def ref(name):
    rows = []
    with open(os.path.join(REF, name)) as f:
        next(f)
        for ln in f:
            if ln.strip():
                a, b = ln.split(',')
                rows.append((float(a), float(b)))
    rows.sort(key=lambda t: t[0])
    r = np.array(rows)
    return r[:, 0], r[:, 1]


for a in sys.argv[1:]:
    vtu, ucsv, pcsv, lab = a.split('=')
    x, uu, pq = centerline(vtu)
    rx_u, ru = ref(ucsv)
    rx_p, rp = ref(pcsv)
    # reference x is already in H units
    ui = np.interp(rx_u, x, uu)
    pi = np.interp(rx_p, x, pq)
    au = np.abs(ui - ru)
    ap = np.abs(pi - rp)
    m2 = rx_u >= 2.0
    m3 = (rx_u >= 3.0) & (rx_u <= 5.0)
    print('== %s   (%s)' % (lab, os.path.basename(vtu)))
    print('   u/U   : L2/max=%.3f%%  L2(x>=2H)/max=%.3f%%  max|du|=%.4f  '
          'mean|du|(2H..)=%.4f' % (100 * np.sqrt((au ** 2).mean()) / ru.max(),
                                   100 * np.sqrt((au[m2] ** 2).mean()) / ru[m2].max(),
                                   au.max(), au[m2].mean()))
    print('           ref  x/H: ' + ' '.join('%6.2f' % v for v in rx_u))
    print('           ours    : ' + ' '.join('%6.3f' % v for v in ui))
    print('           ref     : ' + ' '.join('%6.3f' % v for v in ru))
    print('           diff    : ' + ' '.join('%6.3f' % v for v in (ui - ru)))
    print('   p/rhoU2: L2/max=%.3f%%  max|dp|=%.2f  (ours in 3..5H: %.1f..%.1f, '
          'ref: %.1f..%.1f)' % (100 * np.sqrt((ap ** 2).mean()) / np.abs(rp).max(),
                                ap.max(),
                                pi[(rx_p >= 3) & (rx_p <= 5)].min(),
                                pi[(rx_p >= 3) & (rx_p <= 5)].max(),
                                rp[(rx_p >= 3) & (rx_p <= 5)].min(),
                                rp[(rx_p >= 3) & (rx_p <= 5)].max()))
    print('           ref  x/H: ' + ' '.join('%6.2f' % v for v in rx_p))
    print('           ours    : ' + ' '.join('%6.1f' % v for v in pi))
    print('           ref     : ' + ' '.join('%6.1f' % v for v in rp))
    print('           diff    : ' + ' '.join('%6.1f' % v for v in (pi - rp)))
