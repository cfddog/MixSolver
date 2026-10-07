#!/usr/bin/env python3
"""Compare the coupled channel interface plane (x=100 mm) against the
single-solver full-channel reference (uns_full).

Data sources:
  uns_full/unMesh.vtu     reference, full 0..200 mm channel, pure uns solver
  unMesh_coupled.vtu      coupled run uns side (100..200 mm)
  flow3d.dat              coupled run struct side (0..100 mm, non-dim
                          primitives d,u,v,w,T; parsed with Mesh3d.x dims)

Quantities at the interface:
  * mean u_x and mean p (Pa gauge)
  * u_x(y) profile (20 bins across H=50 mm, z-averaged)
  * struct/uns/reference profile deviation vs the reference
"""
import re
import struct
import sys
import numpy as np

P_REF = 101334.0          # Pa, rho_ref*R*T_ref = 1.177*287.058*300
RHO_REF = 1.177
T_REF = 300.0
U_INF = 34.7224           # m/s, Ma=0.1 * a_ref
H = 0.05
X_IFACE = 0.100
# Both meshes use dx = 2 mm, so the sampling window must be exactly one plane
# wide: [X_IFACE, X_IFACE+DX) picks the single first cell-centre plane
# (x = 101 mm).  A wider window (e.g. 2.5 mm) would also catch x = 103 mm,
# whose pressure is ~220 Pa lower, and manufacture a fictitious jump.
DX = 0.002


def read_vtu_hex(fn):
    """Read an ASCII VTU whose cells are 12 tets per hex (writer order).

    Returns hex centres, pressure (Pa), velocity (m/s), temperature (K).
    """
    txt = open(fn).read()
    ncell = int(re.search(r'NumberOfCells="(\d+)"', txt).group(1))
    pts = np.array(
        re.search(r'Name="Coordinates"[^>]*>(.*?)</DataArray>', txt,
                  re.S).group(1).split(), dtype=float).reshape(-1, 3)
    con = np.array(
        re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt,
                  re.S).group(1).split(), dtype=int).reshape(ncell, -1)
    cc = pts[con[:, 1:]].mean(axis=1)

    def scl(name):
        return np.array(
            re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt,
                      re.S).group(1).split(), dtype=float)

    p = scl('pressure')
    u = scl('velocity').reshape(-1, 3)
    t = scl('temperature')
    nhex = ncell // 12
    assert ncell % 12 == 0, 'not a 12-tet/hex file'
    cx = cc[:, 0].reshape(nhex, 12).mean(axis=1)
    cy = cc[:, 1].reshape(nhex, 12).mean(axis=1)
    cz = cc[:, 2].reshape(nhex, 12).mean(axis=1)
    ph = p.reshape(nhex, 12).mean(axis=1)
    uh = u.reshape(nhex, 12, 3).mean(axis=1)
    th = t.reshape(nhex, 12).mean(axis=1)
    return cx, cy, cz, ph, uh, th


def read_plot3d_dims(fn):
    with open(fn, 'rb') as f:
        rec = f.read(4)
        n = struct.unpack('<i', rec)[0]
        nb = struct.unpack('<i', f.read(4))[0]
        f.read(4)
        f.read(4)
        dims = struct.unpack('<%di' % (3 * nb), f.read(12 * nb))
    return nb, np.array(dims, dtype=int).reshape(nb, 3)


def read_struct_iface(fn, mesh='Mesh3d.x'):
    """Struct flow3d.dat at the i+ interface: returns face u_x(y), p(y),
    z-averaged, SI units (m/s, Pa absolute).

    output_flow writes U(0:nx,0:ny,0:nz) (ghost rows included) where
    bNi = ni (the node count): stored array shape is
    (ni+1)*(nj+1)*(nk+1) in the order i fastest, then j, k, variable.
    Interface i+ (50 cells): inner cell index i=nx-1=ni-1, ghost i=nx=ni.
    """
    nb, dims = read_plot3d_dims(mesh)
    ni, nj, nk = dims[0]          # single block, node counts
    with open(fn, 'rb') as f:
        rec = f.read(4)
        ln = struct.unpack('<i', rec)[0]
        assert ln == 5 * (ni + 1) * (nj + 1) * (nk + 1) * 8, \
               (ln, 5 * (ni + 1) * (nj + 1) * (nk + 1) * 8)
        raw = np.frombuffer(f.read(ln), dtype='<f8')
        f.read(4)
    q = raw.reshape(5, nk + 1, nj + 1, ni + 1)
    d = q[0]
    u = q[1] * U_INF
    T = q[4] * T_REF
    p = d * q[4] * P_REF            # Pa absolute
    # i+ interface: inner cell i=ni-1, ghost i=ni
    uf = 0.5 * (u[:, :, ni - 1] + u[:, :, ni])
    pf = 0.5 * (p[:, :, ni - 1] + p[:, :, ni])
    ui = u[:, :, ni - 1]
    pi = p[:, :, ni - 1]
    # average over interior k layers (1..nk-1) and j cells 1..nj-1
    uf_y = uf[1:nk, 1:nj].mean(axis=0)
    pf_y = pf[1:nk, 1:nj].mean(axis=0)
    ui_y = ui[1:nk, 1:nj].mean(axis=0)
    pi_y = pi[1:nk, 1:nj].mean(axis=0)
    yc = (np.arange(nj - 1) + 0.5) * H / (nj - 1)
    return yc, uf_y, pf_y, ui_y, pi_y


def uns_profile(cx, cy, p, u, x0, x1, ny=20):
    m = (cx >= x0) & (cx < x1)
    ys = np.linspace(0, H, ny + 1)
    yc = 0.5 * (ys[:-1] + ys[1:])
    prof = np.array([
        u[m, 0][(cy[m] >= ya) & (cy[m] < yb)].mean()
        for ya, yb in zip(ys[:-1], ys[1:])])
    pmean = p[m].mean()
    return yc, prof, pmean, m.sum()


def main():
    # ---- reference full channel at x=100 ----
    cx, cy, cz, p, u, t = read_vtu_hex('uns_full/unMesh.vtu')
    yref, uref, pref, ncell = uns_profile(cx, cy, p, u, X_IFACE, X_IFACE + DX)
    print('reference @x=100 (n=%d): mean u=%.3f  mean p=%.1f Pa'
          % (ncell, uref.mean(), pref))

    # ---- coupled uns first layer at its inlet ----
    cx2, cy2, cz2, p2, u2, t2 = read_vtu_hex('unMesh_coupled.vtu')
    xmin = cx2.min()
    yuns, uuns, puns, n2 = uns_profile(cx2, cy2, p2, u2, X_IFACE, X_IFACE + DX)
    print('coupled uns @x=%.3f (n=%d): mean u=%.3f  mean p=%.1f Pa'
          % (xmin, n2, uuns.mean(), puns))

    # ---- coupled struct at i+ interface ----
    ys, uface, pface, uinner, pinner = read_struct_iface('flow3d.dat')
    print('coupled struct face: mean u=%.3f  mean p=%.1f Pa abs (gauge %.1f)'
          % (uface.mean(), pface.mean(), pface.mean() - P_REF))
    print('coupled struct inner: mean u=%.3f  p gauge %.1f'
          % (uinner.mean(), pinner.mean() - P_REF))

    print('\n y(mm)  u_ref  u_struct_face  u_uns_cell0   p_ref  p_struct  p_uns')
    for i in range(len(yref)):
        print('%6.2f %7.2f %8.2f %9.2f   %8.1f %8.1f %8.1f' % (
            yref[i] * 1e3, uref[i], uface[i], uuns[i],
            pref, pface[i] - P_REF, puns))

    # ---- error metrics ----
    du_s = 100 * (uface - uref) / U_INF
    du_u = 100 * (uuns - uref) / U_INF
    print('\nprofile deviation vs reference (%% of U_inf):')
    print('  struct face: max %+.2f%%  RMS %.2f%%  mean-ux %+.2f%%'
          % (np.abs(du_s).max(), np.sqrt((du_s**2).mean()),
             100 * (uface.mean() - uref.mean()) / U_INF))
    print('  uns  cell0 : max %+.2f%%  RMS %.2f%%  mean-ux %+.2f%%'
          % (np.abs(du_u).max(), np.sqrt((du_u**2).mean()),
             100 * (uuns.mean() - uref.mean()) / U_INF))
    print('  pressure: struct %.1f Pa, uns %.1f Pa, ref %.1f Pa (gauge)'
          % (pface.mean() - P_REF, puns, pref))


if __name__ == '__main__':
    main()
