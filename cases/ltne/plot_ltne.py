#!/usr/bin/env python3
"""LTNE 1D validation: numerical T_f / T_s vs the closed-form solution.

Usage: python3 plot_ltne.py <case.vtu> <out.png>

Model (see ltne1d.control header): steady 1D two-phase energy equations
    G*cp*T_f' = h*a*(T_s - T_f),   D*T_s'' = h*a*(T_s - T_f)
with D = (1-eps)*k_s, adiabatic solid at the inlet face and the
conductivity-split wall flux q_s entering the solid at the outlet face.
"""
import re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- case parameters (must match ltne1d.control) --------------------------
eps, ks, kf = 0.3, 16.2, 0.026
G, cp, q   = 2.0, 1005.0, 2.0e6      # mass flux, fluid cp, wall heat flux
ha, L, Tin = 0.3, 50.0, 300.0        # h*a, duct length, inlet temperature
# ---------------------------------------------------------------------------

D  = (1 - eps) * ks
kt = eps * kf + (1 - eps) * ks
qs = q * (1 - eps) * ks / kt         # solid share of the wall flux
b  = ha / (G * cp)
c  = ha / D
disc = np.sqrt(b * b + 4 * c)
lp, lm = (-b + disc) / 2, (-b - disc) / 2
A = (qs / D) / ((1 + lp / b) * (np.exp(lp * L) - np.exp(lm * L)))
B = -A * (lp + b) / (lm + b)

def Tf(x):
    return Tin + A / lp * (np.exp(lp * x) - 1) + B / lm * (np.exp(lm * x) - 1)

def Ts(x):
    return Tf(x) + (A * np.exp(lp * x) + B * np.exp(lm * x)) / b

def read_vtu(fn):
    txt = open(fn).read()
    pts = re.search(r'Name="Coordinates"[^>]*>(.*?)</DataArray>', txt, re.S)
    pts = np.array([float(x) for x in pts.group(1).split()]).reshape(-1, 3)
    con = re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S)
    ncell = int(re.search(r'NumberOfCells="(\d+)"', txt).group(1))
    # records are "4 n1 n2 n3 cen"; the 5th id is the parent cell centroid
    con = np.array([int(x) for x in con.group(1).split()]).reshape(ncell, 5)
    def arr(name, nc):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        return np.array([float(x) for x in m.group(1).split()]).reshape(-1, nc)
    T = arr('temperature', 1).ravel()
    Ts_ = arr('temperature_solid', 1).ravel()
    cen = pts[con[:, 4]]
    return cen, T, Ts_

def main(vtu, png):
    cen, T, Tsg = read_vtu(vtu)
    # bin cell-centred data by x layer (cells of one layer share xc-x)
    xl, Tf_n, Ts_n = np.unique(cen[:, 0]), [], []
    for x in xl:
        m = cen[:, 0] == x
        Tf_n.append(T[m].mean())
        Ts_n.append(Tsg[m].mean())
    Tf_n, Ts_n = np.array(Tf_n), np.array(Ts_n)

    # global energy check: T_f(L)-T_in must equal q/(G*cp)
    dT_ref = q / (G * cp)
    dT_num = Tf_n[-1] - Tin
    print('analytical dT_f(L) = %.2f K   (q/(G*cp) = %.2f K)' %
          (Tf(L) - Tin, dT_ref))
    print('numerical  dT_f(L) = %.2f K   (rel. err %.3f%%)' %
          (dT_num, 100 * (dT_num - dT_ref) / dT_ref))
    err_f = 100 * np.abs(Tf_n - Tf(xl)) / dT_ref
    err_s = 100 * np.abs(Ts_n - Ts(xl)) / (Ts(xl) - Tin).max()
    print('max err T_f (of total rise %.1f K) = %.3f%%   median = %.3f%%'
          % (dT_ref, err_f.max(), np.median(err_f)))
    print('max err T_s (of total rise %.1f K) = %.3f%%   median = %.3f%%'
          % ((Ts(xl) - Tin).max(), err_s.max(), np.median(err_s)))

    xa = np.linspace(0, L, 400)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    ax[0].plot(xa, Tf(xa), 'k-', lw=1.5, label='analytical')
    ax[0].plot(xl, Tf_n, 'o', ms=4, mfc='none', color='tab:red', label='solver')
    ax[0].set_xlabel('x'); ax[0].set_ylabel('$T_f$  [K]')
    ax[0].set_title('fluid temperature')
    ax[1].semilogy(xa, Ts(xa), 'k-', lw=1.5, label='analytical')
    ax[1].semilogy(xl, Ts_n, 'o', ms=4, mfc='none', color='tab:blue', label='solver')
    ax[1].set_xlabel('x'); ax[1].set_ylabel('$T_s$  [K]')
    ax[1].set_title('solid temperature (log scale)')
    ax[2].semilogy(xa, Ts(xa) - Tf(xa), 'k-', lw=1.5, label='analytical')
    ax[2].semilogy(xl, Ts_n - Tf_n, 'o', ms=4, mfc='none',
                   color='tab:green', label='solver')
    ax[2].set_xlabel('x'); ax[2].set_ylabel(r'$\theta=T_s-T_f$  [K]')
    ax[2].set_title('phase temperature difference (log scale)')
    for a in ax:
        a.grid(alpha=0.3); a.legend()
    fig.suptitle(r'LTNE 1D porous slab: air 2 kg/m$^2$s + SS matrix, '
                 r"$q''$=2e6 W/m$^2$ at x=L, $\epsilon$=0.3, $h\cdot a$=0.3")
    fig.tight_layout()
    fig.savefig(png, dpi=130)
    print('wrote', png)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
