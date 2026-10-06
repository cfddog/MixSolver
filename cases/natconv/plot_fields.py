#!/usr/bin/env python3
"""Contour plots (temperature + velocity) from the solver's VTU output.

Usage: python3 plot_fields.py <case.vtu> <out.png> <title>
Cell-centred data are binned onto a uniform 2-D (x,y) grid; the z-extruded
wedge layer makes this an exact 2-D representation.
"""
import re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_vtu(fn):
    txt = open(fn).read()
    pts = re.search(r'Name="Coordinates"[^>]*>(.*?)</DataArray>', txt, re.S)
    pts = np.array([float(x) for x in pts.group(1).split()]).reshape(-1, 3)
    con = re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S)
    ncell = int(re.search(r'NumberOfCells="(\d+)"', txt).group(1))
    # records are "4 n1 n2 n3 cen" (leading node count); 12 tets per hex,
    # the 5th id of every record is the parent hex centroid point.
    con = np.array([int(x) for x in con.group(1).split()]).reshape(ncell, 5)
    def arr(name, nc):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        return np.array([float(x) for x in m.group(1).split()]).reshape(-1, nc)
    vel = arr('velocity', 3)
    T = arr('temperature', 1).ravel()
    cen = pts[con[:, 4]]                 # parent hex centre (exact)
    return cen, vel, T

def bin2d(cen, val, n=200):
    """Bilinear-splat scattered cell data onto a uniform n x n grid."""
    x, y = cen[:, 0], cen[:, 1]
    gx = (np.arange(n) + 0.5) / n
    acc = np.zeros((n, n)); wgt = np.zeros((n, n))
    fx = x * n - 0.5                       # fractional bin coord
    fy = y * n - 0.5
    i0 = np.clip(np.floor(fx).astype(int), 0, n - 2)
    j0 = np.clip(np.floor(fy).astype(int), 0, n - 2)
    ax_ = np.clip(fx - i0, 0, 1)
    ay_ = np.clip(fy - j0, 0, 1)
    for dj in (0, 1):
        for di in (0, 1):
            w = (ay_ if dj else 1 - ay_) * (ax_ if di else 1 - ax_)
            np.add.at(acc, (j0 + dj, i0 + di), val * w)
            np.add.at(wgt, (j0 + dj, i0 + di), w)
    m = wgt > 0
    acc[m] /= wgt[m]
    # fill any still-empty bins from neighbours
    if (~m).any():
        for j, i in zip(*np.where(~m)):
            js, ie = max(j - 1, 0), min(j + 2, n)
            i0_, i1 = max(i - 1, 0), min(i + 2, n)
            nb = acc[js:ie, i0_:i1][wgt[js:ie, i0_:i1] > 0]
            acc[j, i] = nb.mean() if nb.size else 0.0
    return gx, acc

def main(vtu, png, title, nbin=0):
    cen, vel, T = read_vtu(vtu)
    if nbin == 0:                        # bin count ~ sqrt(#cells)
        nbin = max(20, int(round(np.sqrt(len(T)))))
    xs, Tg = bin2d(cen, T - 300.5, n=nbin)
    sp = np.hypot(vel[:, 0], vel[:, 1])
    ns = max(20, nbin // 3)
    xs60, Ug = bin2d(cen, vel[:, 0], n=ns)
    _, Vg = bin2d(cen, vel[:, 1], n=ns)
    _, Sg = bin2d(cen, sp, n=nbin)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))
    cf = ax[0].contourf(xs, xs, Tg, levels=40, cmap='RdBu_r')
    ax[0].set_title('Temperature  $T-T_{ref}$  [K]')
    ax[0].set_aspect('equal')
    plt.colorbar(cf, ax=ax[0], shrink=0.9)
    cs = ax[1].contourf(xs, xs, Sg, levels=40, cmap='viridis')
    ax[1].streamplot(xs60, xs60, Ug, Vg, color='w', density=0.9,
                     linewidth=0.7, arrowsize=0.9)
    ax[1].set_title('Velocity magnitude  |u|  + streamlines')
    ax[1].set_aspect('equal')
    plt.colorbar(cs, ax=ax[1], shrink=0.9)
    for a in ax:
        a.set_xlabel('x'); a.set_ylabel('y')
    fig.suptitle(title, y=0.99)
    fig.tight_layout()
    fig.savefig(png, dpi=150)
    print('wrote', png)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3],
         int(sys.argv[4]) if len(sys.argv) > 4 else 0)
