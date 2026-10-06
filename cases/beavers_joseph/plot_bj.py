#!/usr/bin/env python3
# Beavers-Joseph validation against the EXACT two-layer reference.
#
# The discrete BJ flux converges to a "stress continuity + velocity jump"
# interface model, NOT to the classical BJ slip formula (which assumes a
# pure-Darcy porous side).  With porosity=1 the Brinkman term is active in
# the porous layer, so the continuum reference is the coupled ODE system
# (y measured from the interface, fluid y in [0,H], porous y in [-H,0]):
#
#   fluid : mu u_f'' + fb = 0,                  u_f(H)  = 0
#   porous: mu u_p'' - (mu/K) u_p + fb = 0,     u_p(-H) = 0
#   interface: mu u_f'(0) = mu u_p'(0) = tau
#     alpha>0 : tau = mu (alpha/lam) (u_f(0) - u_p(0))   (stress-jump)
#     alpha=0 : u_f(0) = u_p(0)                          (continuity)
#
# with lam = sqrt(K), u_D = K fb/mu:
#   u_f = -fb y^2/(2 mu) + a1 y + a2
#   u_p = u_D + E exp(y/lam) + F exp(-y/lam)
# -> 4x4 linear solve for [a1, a2, E, F].
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LX, LY, NX, NY = 0.10, 0.010, 100, 40
H = LY / 2
MU, K, FB = 1.846e-5, 1.0e-6, 0.05
LAM = np.sqrt(K)
UD = K * FB / MU
SIGMA = H / LAM


def load(fn):
    s = open(fn).read()

    def arr(name, text):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii">(.*?)</DataArray>'
                      % re.escape(name), text, re.S)
        return np.fromstring(m.group(1), sep=' ')

    sec = re.search(r'<Cells>(.*?)</Cells>', s, re.S).group(1)
    conn = arr('connectivity', sec).astype(int)
    off = arr('offsets', sec).astype(int)
    cd_ = re.search(r'<CellData>(.*?)</CellData>', s, re.S).group(1)
    uvw = arr('velocity', cd_).reshape(-1, 3)
    pts = arr('Coordinates', re.search(r'<Points>(.*?)</Points>', s, re.S)
              .group(1)).reshape(-1, 3)
    st = np.concatenate(([0], off[:-1]))
    xyc = np.zeros((len(off), 2))
    for i in range(len(off)):
        ids = conn[st[i]:off[i]]
        if ids[0] == 4 and len(ids) == 5:
            ids = ids[1:]
        xyc[i] = pts[ids].mean(0)[:2]
    ix = np.round((xyc[:, 0] - LX / NX / 2) / (LX / NX)).astype(int)
    iy = np.round((xyc[:, 1] - LY / NY / 2) / (LY / NY)).astype(int)
    return ix, iy, uvw


def reference(alpha):
    # unknowns x = [a1, a2, E, F]
    em, ep = np.exp(-H / LAM), np.exp(H / LAM)
    A = np.zeros((4, 4))
    b = np.zeros(4)
    # u_f(H) = 0
    A[0] = [H, 1.0, 0.0, 0.0]
    b[0] = FB * H**2 / (2 * MU)
    # u_p(-H) = 0
    A[1] = [0.0, 0.0, em, ep]
    b[1] = -UD
    # stress continuity: a1 = (E - F)/lam
    A[2] = [1.0, 0.0, -1.0 / LAM, 1.0 / LAM]
    # interface velocity condition
    if alpha > 0.0:
        # a1 = (alpha/lam)(a2 - uD - E - F)
        A[3] = [1.0, -alpha / LAM, alpha / LAM, alpha / LAM]
        b[3] = -alpha / LAM * UD
    else:
        # a2 = uD + E + F
        A[3] = [0.0, 1.0, -1.0, -1.0]
        b[3] = UD
    a1, a2, E, F = np.linalg.solve(A, b)

    def uf(y):
        return -FB * y**2 / (2 * MU) + a1 * y + a2

    def up(y):
        return UD + E * np.exp(y / LAM) + F * np.exp(-y / LAM)

    return uf, up, a2  # a2 = u_f(0) interface speed (fluid side)


def classical_uB(alpha):
    # classical BJ slip speed (pure-Darcy porous side), for comparison only
    return (FB * H**2 / (2 * MU) + alpha * SIGMA * UD) / (1 + alpha * SIGMA)


runs = [('bj_a0.vtu', 0.0, 'C0', 'continuity (alpha=0)'),
        ('bj_a1.vtu', 1.0, 'C1', 'BJ alpha=1'),
        ('bj_a2.vtu', 2.0, 'C2', 'BJ alpha=2')]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
yc = (np.arange(NY) + 0.5) * LY / NY          # cell centres from bottom wall
yf = np.linspace(0, H, 120)
yp = np.linspace(-H, 0, 120)

print('%6s %11s %11s %11s %11s %11s %9s'
      % ('alpha', 'uF1_num', 'uF1_ref', 'uD_num', 'uD_ref', 'uB_class',
         'rms_err'))
dy = LY / NY
# porous-bed interior rows: >= lam away from the bottom wall and the interface
r0, r1 = int(np.ceil(LAM / dy)), NY // 2 - int(np.ceil(LAM / dy))
for fn, alpha, col, lab in runs:
    ix, iy, uvw = load(fn)
    colx = NX // 2
    ucol = np.array([uvw[(ix == colx) & (iy == j), 0].mean()
                     for j in range(NY)])
    uf, up, uif_ref = reference(alpha)
    # reference at cell centres (y relative to the interface)
    yrel = yc - H
    uref = np.where(yrel >= 0.0, uf(np.maximum(yrel, 0.0)),
                    up(np.minimum(yrel, 0.0)))
    rms = np.sqrt(np.mean((ucol - uref)**2)) / max(uref) * 100
    uD_num = ucol[r0:r1].mean()               # porous plug away from layers
    uF1_num = ucol[NY // 2]                   # first fluid cell (y=dy/2)
    uF1_ref = uf(dy / 2)
    print('%6.1f %11.4e %11.4e %11.4e %11.4e %11.4e %8.2f%%'
          % (alpha, uF1_num, uF1_ref, uD_num, UD, classical_uB(alpha), rms))

    axes[0].plot(ucol, yc, 'o', ms=3, color=col, label=lab)
    axes[0].plot(up(yp), yp + H, '-', color=col, lw=1.2)
    axes[0].plot(uf(yf), yf + H, '-', color=col, lw=1.2)

    axes[1].plot(ucol, yrel, 'o', ms=4, color=col, label=lab)
    axes[1].plot(up(yp), yp, '-', color=col, lw=1.2)
    axes[1].plot(uf(yf), yf, '-', color=col, lw=1.2)
    if alpha > 0.0:
        axes[1].axvline(classical_uB(alpha), color=col, ls=':', lw=1,
                        label='classical BJ a=%g' % alpha)

axes[0].axhline(H, color='k', ls='--', lw=0.8)
axes[0].set_xlabel('u (m/s)')
axes[0].set_ylabel('y (m)')
axes[0].set_title('Full channel profile at x=L/2')
axes[0].legend(fontsize=8)
axes[0].grid(alpha=0.3)

axes[1].axhline(0.0, color='k', ls='--', lw=0.8)
axes[1].set_xlabel('u (m/s)')
axes[1].set_ylabel('y - H (m)')
axes[1].set_title('Interface zoom (lines: two-layer reference)')
axes[1].legend(fontsize=7)
axes[1].grid(alpha=0.3)
axes[1].set_ylim(-0.002, 0.005)

fig.tight_layout()
fig.savefig('images/beavers_joseph.png', dpi=140)
print('wrote images/beavers_joseph.png')
