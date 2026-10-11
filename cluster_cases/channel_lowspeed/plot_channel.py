#!/usr/bin/env python3
"""Planar Poiseuille channel validation: profile, entrance length, pressure.

Usage: python3 plot_channel.py <case.vtu> <out.png> [Re_H]

The mesh is 300 (x) x 100 (y, wall-clustered) x 1 (z) hex cells.  VTU hex
cells are decomposed into 12 tets with records "4 n1 n2 n3 cen": the 5th
id is the parent hex centroid, and cell data are replicated 12 times, so
unique centroid ids recover one record per solver cell.

Reference (steady, incompressible, laminar, full height H = y_top - y_bot):
    u(y) = 6 U_bulk eta(1-eta), eta = y/H ; u_max/U = 1.5
    dp/dx = -12 mu U / H^2 ; L_e/H ~ 0.05 Re_H
"""
import re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- case parameters (match channel_Re*.control) ---
rho, mu = 1.177, 1.846e-5
H = 10.0
Re = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0
U = Re * mu / (rho * H)
dpdx_an = -12.0 * mu * U / H**2

def read_vtu(fn):
    txt = open(fn).read()
    pts = re.search(r'Name="Coordinates"[^>]*>(.*?)</DataArray>', txt, re.S)
    pts = np.array([float(x) for x in pts.group(1).split()]).reshape(-1, 3)
    ncell = int(re.search(r'NumberOfCells="(\d+)"', txt).group(1))
    con = re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S)
    con = np.array([int(x) for x in con.group(1).split()]).reshape(ncell, 5)

    def vec(name):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        return np.array([float(x) for x in m.group(1).split()]).reshape(-1, 3)

    def scl(name):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        return np.array([float(x) for x in m.group(1).split()]).reshape(-1, 1).ravel()

    vel, pr = vec('velocity'), scl('pressure')
    cenid, inv = np.unique(con[:, 4], return_inverse=True)
    cen = pts[cenid]
    # replicated values are identical across the 12 tets; take first occurrence
    u = np.zeros(len(cenid)); v = np.zeros(len(cenid)); w = np.zeros(len(cenid))
    p = np.zeros(len(cenid))
    seen = np.zeros(len(cenid), bool)
    for k in range(len(inv)):
        j = inv[k]
        if not seen[j]:
            u[j], v[j], w[j] = vel[k]
            p[j] = pr[k]
            seen[j] = True
    return cen, u, v, w, p

def layer_weights(y):
    """trapezoidal integration weights for clustered cell-center layers."""
    order = np.argsort(y)
    ys = y[order]
    edges = np.empty(len(ys) + 1)
    edges[1:-1] = 0.5 * (ys[:-1] + ys[1:])
    edges[0] = max(ys[0] - (ys[1] - ys[0]) / 2, 0.0)
    edges[-1] = min(ys[-1] + (ys[-1] - ys[-2]) / 2, H)
    wt = np.diff(edges)
    w = np.zeros_like(y)
    w[order] = wt
    return w

def main(vtu, png):
    cen, u, v, w, p = read_vtu(vtu)
    X, Y = cen[:, 0], cen[:, 1]
    xlvl = np.unique(np.round(X, 9))

    # ---- per-x-layer bulk averages ----
    Ubulk = []; Ucl = []; play = []
    wt0 = None
    for x in xlvl:
        m = np.abs(X - x) < 1e-6
        y_, u_ = Y[m], u[m]
        wt = layer_weights(y_)
        Ubulk.append(np.sum(u_ * wt) / H)
        play.append(np.sum(p[m] * wt) / H)
        Ucl.append(u_[np.argmin(np.abs(y_ - H / 2))])
    Ubulk, Ucl, play = map(np.array, (Ubulk, Ucl, play))

    # ---- mass conservation across stations ----
    print('U_bulk target = %.6e  (Re_H=%.0f)' % (U, Re))
    print('U at x=0.5 / 50.5 / 299.5 = %.6e %.6e %.6e' %
          (Ubulk[0], Ubulk[50], Ubulk[-1]))
    print('mass drift in->out = %.3f%%' % (100 * (Ubulk[-1] - Ubulk[0]) / Ubulk[0]))

    # ---- developed-profile comparison (last 30 length used for stats) ----
    eta = np.linspace(0, 1, 200)
    ua = 6 * U * eta * (1 - eta)
    stations = [5.5, 25.5, 50.5, 149.5, 299.5]
    prof = {}
    for xs in stations:
        k = int(np.argmin(np.abs(xlvl - xs)))
        m = np.abs(X - xlvl[k]) < 1e-6
        y_, u_ = Y[m], u[m]
        o = np.argsort(y_)
        prof[xs] = (y_[o], u_[o])

    yd, ud = prof[299.5]
    uan = 6 * U * (yd / H) * (1 - yd / H)
    err = np.abs(ud - uan) / U
    print('developed profile (x=299.5): max err = %.3f%% of U, L1 err = %.3f%%'
          % (100 * err.max(), 100 * err.mean()))
    print('u_max/U at outlet = %.4f  (analytic 1.5)' % (ud.max() / U))

    # ---- entrance length: |Ucl - 1.5U| within 1% ----
    band = np.abs(Ucl - 1.5 * U) / U < 0.01
    Le = xlvl[np.argmax(band)] if band.any() else np.nan
    print('entrance length Le = %.1f (%.2f H ; estimate 0.05ReH = %.1f)'
          % (Le, Le / H, 0.05 * Re * H / H))

    # ---- developed pressure gradient (linear fit x in [200,295]) ----
    mf = (xlvl >= 200) & (xlvl <= 295)
    slope, intercept = np.polyfit(xlvl[mf], play[mf], 1)
    print('dp/dx fitted [200,295] = %.6e   analytic = %.6e   err %.2f%%'
          % (slope, dpdx_an, 100 * (slope - dpdx_an) / dpdx_an))
    print('total pressure rise in->out = %.4e' % (play[0] - play[-1]))

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
    for xs in stations:
        yy, uu = prof[xs]
        ax[0].plot(uu / U, yy / H, 'o', ms=3, mfc='none',
                   label='x=%g' % xs)
    ax[0].plot(ua / U, eta, 'k-', lw=1.4, label='analytic')
    ax[0].set_xlabel(r'$u/U_{bulk}$'); ax[0].set_ylabel(r'$y/H$')
    ax[0].set_title('streamwise profiles')

    ax[1].plot(xlvl / H, Ucl / U, '.', ms=3)
    ax[1].axhline(1.5, color='k', ls='--', lw=1)
    ax[1].axvline(Le / H, color='r', ls=':', lw=1)
    ax[1].set_xlabel(r'$x/H$'); ax[1].set_ylabel(r'$u_{cl}/U_{bulk}$')
    ax[1].set_title('centerline development, $L_e$=%.1f (%.2f H)' % (Le, Le / H))

    ax[2].plot(xlvl, play, '.', ms=3, label='solver (layer mean)')
    ax[2].plot(xlvl[mf], slope * xlvl[mf] + intercept, 'r-', lw=1.2,
               label='fit [200,295]\nslope err %.2f%%' %
               (100 * (slope - dpdx_an) / dpdx_an))
    ax[2].set_xlabel('$x$'); ax[2].set_ylabel('$p$')
    ax[2].set_title(r'pressure, analytic $dp/dx$=%.3e' % dpdx_an)
    for a in ax:
        a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.suptitle(r'Planar Poiseuille channel, $Re_H$=%.0f, air 300 K, '
                 r'mass-flow-inlet + pressure-outlet' % Re)
    fig.tight_layout()
    fig.savefig(png, dpi=130)
    print('wrote', png)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
