#!/usr/bin/env python3
"""Cylinder Re=40 scene check from a bin/mixnsolver output .vtu.

Usage:  python3 check_cylinder.py <case.vtu> [reference_cylinder_cp.dat]

Case: 2D laminar cylinder, D = 1, U_inf = (1,0,0), Re = rho*U*D/mu = 40,
far-field at r = 10, strictly 2D (one cell in z -> per-unit-span forces).

What is printed
  * cell / node counts, max|u| (must be ~U_inf), p range
  * Cp = p / (0.5 rho U^2) on the innermost radial cell ring (wall-adjacent):
    Cp_min (front stagnation, expect ~ +1), Cp_max (base suction)
  * separation angle measured from the front stagnation point (expect ~53 deg)
  * recirculation length on the wake axis in diameters (expect ~2.3 D)
  * pressure-drag coefficient CD_p (expect ~1.0; the literature TOTAL CD at
    Re=40 is ~1.5, the difference is the friction part, not resolved here)
  * optional: |dCp| statistics against the archived Cp(theta) baseline

VTU layout of this code: points = mesh nodes then one cell-centroid point per
cell; cells = tetrahedra from a face fan + the cell centroid, written as
"n_pts  i1 i2 i3  centroid_id" (5 integers per tet, see vtk_write in
src/unstructured/mod_uns_output.f90).  Cell data are replicated on every
sub-tet, so unique centroid ids recover one value per solver cell.

All metrics are regression indicators (reproduce them on the cluster); the
wall pressure is cell-centred, i.e. it carries an O(dy^2) offset on this
coarse 80 x 49 mesh.
"""
import re
import sys

import numpy as np

R_CYL = 0.5          # cylinder radius (D = 1)
U_INF = 1.0
RHO = 1.0
Q = 0.5 * RHO * U_INF**2


def read_vtu(fn):
    txt = open(fn).read()

    def arr(name, cols=1, dtype=float):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % re.escape(name),
                      txt, re.S)
        if m is None:
            raise SystemExit('check_cylinder: no "%s" array in %s' % (name, fn))
        a = np.array([dtype(x) for x in m.group(1).split()])
        return a.reshape(-1, cols) if cols > 1 else a

    pts = arr('Coordinates', 3)
    n_tets = int(re.search(r'NumberOfCells="(\d+)"', txt).group(1))
    con = arr('connectivity', 1, int)          # flat: node ids of every tet
    con = con.reshape(n_tets, -1)
    if con.shape[1] == 5 and con[:, 0].max() == 4:   # old files: legacy VTK
        con = con[:, 1:]                             # leading vertex-count token
    if con.shape != (n_tets, 4):
        raise SystemExit('check_cylinder: unexpected connectivity layout')
    pres = arr('pressure')
    vel = arr('velocity', 3)

    # one value per solver cell (first sub-tet of each centroid group)
    cen_id = con[:, -1]
    first = np.concatenate(([True], cen_id[1:] != cen_id[:-1]))
    idx = np.nonzero(first)[0]
    base = cen_id.min()                        # 0-based id of cell-1 centroid
    xyz = pts[cen_id[idx], :]                  # cell centroids
    return (xyz, vel[idx, :], pres[idx], cen_id[idx] - base + 1, n_tets)


def radial_ring(xc, yc, ux, p):
    """Innermost-wall ring: one cell per angular bin (80 bins), min radius."""
    th = np.degrees(np.arctan2(yc, xc))
    rr = np.hypot(xc, yc)
    bin_ = np.clip(((th + 180.0) / 360.0 * 80.0).astype(int), 0, 79)
    sel = []
    for b in range(80):
        m = bin_ == b
        if m.any():
            j = np.nonzero(m)[0][np.argmin(rr[m])]
            sel.append(j)
    sel = np.array(sel, dtype=int)
    o = np.argsort(th[sel])
    return th[sel][o], ux[sel][o], p[sel][o], rr[sel][o]


def wake_axis(xc, yc, ux):
    """Cell nearest to y = 0 in every x bin -> centreline u_x(x)."""
    xs, us = [], []
    for x0 in np.arange(0.5, 6.0, 0.05):
        m = (xc > x0 - 0.025) & (xc < x0 + 0.025)
        if m.any():
            j = np.nonzero(m)[0][np.argmin(np.abs(yc[m]))]
            xs.append(xc[j].mean())
            us.append(ux[j].mean())
    return np.array(xs), np.array(us)



def main(fn, ref_fn=None):
    xyz, vel, pres, cid, n_tets = read_vtu(fn)
    xc, yc = xyz[:, 0], xyz[:, 1]
    ux = vel[:, 0]
    rr = np.hypot(xc, yc)

    print('file            : %s' % fn)
    print('solver cells    : %d    tets: %d' % (cid.max(), n_tets))
    print('cell r range    : %.4f .. %.4f  (cylinder R=%.2f, far field r=10)'
          % (rr.min(), rr.max(), R_CYL))
    print('x range         : %.3f .. %.3f' % (xc.min(), xc.max()))
    print('max|u| / U_inf  : %.6f    (ideal 1.0)'
          % (np.hypot(vel[:, 0], vel[:, 1]).max() / U_INF))
    print('p range         : %.6e .. %.6e     Cp range %.4f .. %.4f'
          % (pres.min(), pres.max(), 2 * pres.min() / Q, 2 * pres.max() / Q))

    th_w, ux_w, p_w, r_w = radial_ring(xc, yc, ux, pres)
    cp_w = p_w / Q
    # phi = angle measured FROM the front stagnation point (far field is +x,
    # so the front stagnation is at theta = 180 deg, the wake axis at phi=180)
    phi = (th_w % 360.0) - 180.0
    i_f = int(np.argmin(np.abs(phi)))
    i_r = int(np.argmax(np.abs(phi)))
    print('wall ring cells : %d   (mean r = %.4f, dr_first = %.4f)'
          % (len(th_w), r_w.mean(), r_w.min() - R_CYL))
    print('Cp front (phi=0): %.4f   (expect ~ +1.5, potential-flow limit 1.0)'
          % cp_w[i_f])
    print('Cp base  (phi=180): %.4f  (suction in the dead-water zone)'
          % cp_w[i_r])
    print('Cp suction peak : %.4f at phi = %+.1f deg  (expect ~ -0.97)'
          % (cp_w.min(), phi[int(np.argmin(cp_w))]))

    # separation: near-wall u_x sign change, reported as |phi| (from the front
    # stagnation point) and equivalently as (180 - |phi|) from the wake axis
    o = np.argsort(phi)
    ph, uw = phi[o], ux_w[o]
    sgn = np.sign(uw)
    chg = np.nonzero(sgn[1:] * sgn[:-1] < 0)[0]
    if len(chg):
        for i in chg:
            print('separation      : |phi| = %.1f deg from front stagnation'
                  '  (%.1f deg from wake axis)'
                  % (abs(ph[i]), 180.0 - abs(ph[i])))
    else:
        print('separation      : no near-wall reversal found (attached)')

    # recirculation length on the wake axis
    xs, us = wake_axis(xc, yc, ux)
    neg = us < 0.0
    if neg.any():
        print('wake recirc.    : L_w = %.4f  (%.3f D)   max reverse u_x = %.5f'
              % (xs[neg].max(), xs[neg].max() / (2 * R_CYL), us.min()))
    else:
        print('wake recirc.    : none (attached wake)')

    # pressure drag from the wall ring:  CD = -sum(Cp n_x) R dtheta / D
    dth = np.diff(np.concatenate((th_w, th_w[:1] + 360.0)))
    cd_p = -np.sum(cp_w * np.cos(np.radians(th_w)) * R_CYL
                   * np.radians(dth)) / (2 * R_CYL)
    print('CD_pressure     : %.4f   (literature total CD(Re=40) ~ 1.5)'
          % cd_p)

    if ref_fn:
        ref = np.loadtxt(ref_fn, comments='#')
        rth, rcp = ref[:, 0], ref[:, 1]
        rth = np.where(rth > 180.0, rth - 360.0, rth)
        o = np.argsort(rth)
        rth, rcp = rth[o], rcp[o]
        thx = np.concatenate((rth - 360.0, rth, rth + 360.0))
        cpx = np.tile(rcp, 3)
        d = cp_w - np.interp(th_w, thx, cpx)
        print('vs %s (%d pts): mean|dCp| = %.4f   max|dCp| = %.4f'
              % (ref_fn.split('/')[-1], len(rth), np.abs(d).mean(),
                 np.abs(d).max()))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
