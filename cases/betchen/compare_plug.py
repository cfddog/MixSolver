#!/usr/bin/env python3
# Betchen PLUG comparison: centerline velocity u/U and pressure p/(rho U^2)
# along x/H against the provided reference (plug_1 = Da=1e-2, plug_2 =
# Da=1e-3).  Reference inlet is parabolic fully-developed u=6U y/H(1-y/H),
# whereas the current solver's velocity-inlet is zone-uniform -- so the
# upstream fluid section will differ (ref centerline 1.5U vs ours ~1.0U);
# the porous region and downstream are the meaningful comparison.
import re, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H = 0.01
RHO, MU = 1.177, 1.846e-5
U0 = 1.5684e-3                       # our uniform inlet (mean) velocity
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


def centerline(cc, ucn, pp, ny=21):
    # cells on the centerline row (the row containing y=H/2)
    j = np.abs(cc[:, 1] - H / 2) < (H / ny)
    xs = cc[j, 0] / H
    uu = ucn[j] / U0
    pq = pp[j] / (RHO * U0 ** 2)
    o = np.argsort(xs)
    return xs[o], uu[o], pq[o]


def get_ref(name):
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


cases = [
    ('plug_dae2.vtu', 'plug_1_u.csv', 'plug_1_p.csv', 'PLUG Da=1e-2 Re=1', 'C0'),
    ('plug_dae3.vtu', 'plug_2_u.csv', 'plug_2_p.csv', 'PLUG Da=1e-3 Re=1', 'C3'),
]

fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
for vtu, ucsv, pcsv, lab, col in cases:
    cc, ucn, pp = load(vtu)
    x, uu, pq = centerline(cc, ucn, pp)
    rx_u, ru = get_ref(os.path.join(REF, ucsv))
    rx_p, rp = get_ref(os.path.join(REF, pcsv))
    # interpolate mine onto ref x
    ui = np.interp(rx_u, x, uu)
    pi = np.interp(rx_p, x, pq)
    l2u = np.sqrt(np.mean((ui - ru) ** 2)) / (np.abs(ru).max())
    l2p = np.sqrt(np.mean((pi - rp) ** 2)) / (np.abs(rp).max() + 1e-30)
    # exclude the developing inlet region for the velocity metric (x/H>=2)
    m = rx_u >= 2.0
    l2u2 = np.sqrt(np.mean((ui[m] - ru[m]) ** 2)) / (np.abs(ru[m]).max())
    print('== %s' % lab)
    print('   velocity u/U: L2(mx)/maxU=%6.3f%%   L2(x/H>=2)/maxU=%6.3f%%   '
          'ref-full=%6.3f%%' % (100*l2u, 100*l2u2, 100*(ui.max()-ru.min())))
    print('   pressure  p/(rhoU^2): L2/max|prev|=%6.3f%%   '
          'my range[%.1f, %.1f] ref range[%.1f, %.1f]'
          % (100*l2p, pi.min(), pi.max(), rp.min(), rp.max()))
    # solver = solid circles with no markers edges; ref = heavy dashed crosses,
    # same hue per case so the pair is visually grouped but the two roles differ
    axs[0].plot(x, uu, '-o', ms=5, lw=1.8, color=col, fillstyle='none',
                label='SIMPLE ' + lab)
    axs[0].plot(rx_u, ru, '--.', lw=2.2, ms=11, color=col, label='ref ' + lab)
    axs[1].plot(x, pq, '-o', ms=5, lw=1.8, color=col, fillstyle='none',
                label='SIMPLE ' + lab)
    axs[1].plot(rx_p, rp, '--.', lw=2.2, ms=11, color=col, label='ref ' + lab)

axs[0].set_xlabel('x / H'); axs[0].set_ylabel('u/U (centerline)')
axs[0].set_title('Centerline velocity along flow')
axs[0].grid(alpha=0.3); axs[0].legend(fontsize=7, ncol=2)
axs[1].set_xlabel('x / H'); axs[1].set_ylabel('p / (rho U^2)')
axs[1].set_title('Centerline pressure along flow')
axs[1].grid(alpha=0.3); axs[1].legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig('images/plug_compare_ref.png', dpi=150)
print('wrote images/plug_compare_ref.png')