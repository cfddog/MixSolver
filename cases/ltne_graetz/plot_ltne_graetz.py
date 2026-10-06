#!/usr/bin/env python3
# LTNE Graetz validation: fluid/solid column profiles vs the coupled
# two-equation Graetz quadratic eigenproblem (semi-analytic reference).
#
# Steady plug flow (superficial U) in a porous half-channel, eta=y/a in
# [0,1] (adiabatic centreline eta=0, fixed-T wall eta=1), xi=x/a.  Both
# phases retain axial conduction (epsilon*kf for the fluid, (1-eps)*ks
# for the solid), matching the solver assembly:
#   (1/Pex) theta_f,xixi + (1/Pe) L22 theta_f
#       = theta_f,xi + Bi (theta_f-theta_s)
#   Lams (theta_s,xixi + L22 theta_s) = Bi (theta_s - theta_f)
# 1/Pe = Kfy/(rho cp U a) transverse, 1/Pex = Kfx/(...) axial,
# Bi = H a/(rho cp U), Lams = Ds/(rho cp U a),
# Kfy = eps*kf + rho cp disp_t U, Kfx = eps*kf + rho cp disp_l U,
# Ds = (1-eps) k_s, H = h_sf a_sf.
#
# The quadratic eigenproblem is written as a first-order block system
# y'=M y, y=[f,s,chi_f,chi_s]; decaying modes exp(-lam xi) are retained
# and the entrance is fixed by theta_f(0)=1 and chi_s(0)=0 (the solver
# enforces Tf=Tin at the inlet and treats the solid inlet as adiabatic).
# theta = (T - Tw)/(T0 - Tw); x* = x/(Pe a^2) = kf x/(U a^2),
# modes decay as exp(-lam Pe x*).
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LX, LY, NX, NY = 0.45, 0.005, 180, 40
A_H = LY
RHO, CP, KF = 1.177, 1005.0, 0.026
EPS, KS, U, H = 0.4, 1.0, 1.0, 2.0e5
T0, TW = 300.0, 310.0
DS = (1 - EPS) * KS


def load(fn, nx=NX, ny=NY):
    s = open(fn).read()

    def arr(name, text):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii">(.*?)</DataArray>'
                      % re.escape(name), text, re.S)
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
    xyc = np.zeros((len(off), 2))
    for i in range(len(off)):                 # hex exported as 12 tets
        ids = conn[st[i]:off[i]]
        if ids[0] == 4 and len(ids) == 5:
            ids = ids[1:]
        xyc[i] = pts[ids].mean(0)[:2]
    ix = np.round((xyc[:, 0] - LX / nx / 2) / (LX / nx)).astype(int)
    iy = np.round((xyc[:, 1] - LY / ny / 2) / (LY / ny)).astype(int)
    return ix, iy, Tf, Ts


# Both phases keep their axial conduction (epsilon*kf for the fluid,
# (1-eps)*ks for the solid).  In (xi,eta), theta=(T-Tw)/(T0-Tw):
#   (1/Pe)(theta_f,xi xi + L theta_f) = theta_f,xi + Bi(theta_f - theta_s)
#   Lams (theta_s,xi xi + L theta_s) = Bi(theta_s - theta_f)
# Modes exp(mu xi) form a quadratic eigenproblem; write it as the
# first-order block system  y' = M y,  y=[f, s, chi_f, chi_s],
# chi= d()/dxi, keep the 2m decaying modes (Re mu < 0, lam=-mu) and
# fix the entrance with  theta_f(0)=1, chi_s(0)=0 (solid face adiabatic;
# exactly the code's inlet treatment).

def eigenmodes(disp_t, ngrid=180, nmodes=200, disp_l=0.0):
    # transverse fluid conductance: epsilon*kf + rho cp disp_t U;
    # axial: epsilon*kf + rho cp disp_l U (disp_l=0 here)
    Kfy = EPS * KF + RHO * CP * disp_t * U
    Kfx = EPS * KF + RHO * CP * disp_l * U
    kappa = Kfy / (RHO * CP)
    Pe = U * A_H / kappa                 # transverse Peclet (x* scaling)
    Pex = U * A_H / (Kfx / (RHO * CP))  # axial-conduction Peclet
    Bi = H * A_H / (RHO * CP * U)
    Lams = DS / (RHO * CP * U * A_H)

    m = ngrid
    h = 1.0 / m
    eta = (np.arange(m) + 0.5) * h           # cell-centre-like reference grid
    L2 = np.zeros((m, m))
    for i in range(1, m - 1):
        L2[i, i - 1:i + 2] = [1, -2, 1]
    L2[0] = np.array([-6, 7, -1] + [0] * (m - 3)) / 3.0   # Neumann eta=0
    L2[m - 1, m - 2], L2[m - 1, m - 1] = 1, -2            # phi(1)=0
    L2 /= h**2

    I = np.eye(m)
    Z = np.zeros((m, m))
    M = np.block([
        [Z, Z, I, Z],
        [Z, Z, Z, I],
        [Pex * Bi * I - (Pex / Pe) * L2, -Pex * Bi * I, Pex * I, Z],
        [-(Bi / Lams) * I, -L2 + (Bi / Lams) * I, Z, Z]])
    w, V = np.linalg.eig(M)
    dec = w.real < -1e-8
    w, V = w[dec], V[:, dec]
    lam = -w.real
    Vf = V[0:m, :].real
    Vs = V[m:2 * m, :].real
    # scale modes by their largest phase amplitude for conditioning
    sc = np.maximum(np.abs(Vf).max(0), np.abs(Vs).max(0))
    Vf, Vs = Vf / sc, Vs / sc
    order = np.argsort(lam)
    lam, Vf, Vs = lam[order][:nmodes], Vf[:, order][:, :nmodes], \
        Vs[:, order][:, :nmodes]
    # entrance coefficients: theta_f(0)=1, chi_s(0)=sum c mu phi_s=0
    chi_s = (-lam) * Vs
    A0 = np.vstack([Vf, chi_s])
    b0 = np.concatenate([np.ones(m), np.zeros(m)])
    c, *_ = np.linalg.lstsq(A0, b0, rcond=None)
    return kappa, Pe, Bi, Lams, eta, lam, c, Vf, Vs


def profiles(disp_t):
    return eigenmodes(disp_t)


def theta_at(xs_star, ref):
    kappa, Pe, eta, lam, c, Vf, Vs = ref[0], ref[1], ref[4], ref[5], \
        ref[6], ref[7], ref[8]
    ef = np.exp(-lam * Pe * xs_star)
    return Vf @ (c * ef), Vs @ (c * ef)


runs = [('ltne_grz_d.vtu', 2.5e-4, 'C3', 'LTNE + disp_t=2.5e-4 m'),
        ('ltne_grz0.vtu', 0.0, 'C0', 'LTNE, disp_t=0')]
# common physical stations (m); per-case x* and columns are derived below
xphys = [0.02, 0.05, 0.10, 0.20, 0.35]
pcols = ['C0', 'C1', 'C2', 'C3', 'C4']

fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
eta_d = (np.arange(NY) + 0.5) / NY
for row, (fn, al, base, title) in enumerate([
        ('ltne_grz_d.vtu', 2.5e-4, 'C3', 'with transverse dispersion'),
        ('ltne_grz0.vtu', 0.0, 'C0', 'conduction only (disp_t=0)')]):
    ix, iy, Tf, Ts = load(fn)
    ixx, iyx, Tfx, Tsx = load(fn.replace('.vtu', '_xref.vtu'), 360, 40)
    ref = profiles(al)
    kappa, Pe, Bi, Lams, eta, lam, c, Vf, Vs = ref
    print('%-26s Pe=%7.2f Bi=%.3f Lams=%.4f nmodes=%d'
          % (title, Pe, Bi, Lams, len(lam)))
    # entrance check: solid profile at the first cell column vs the
    # eigen expansion with chi_s(0)=0 (adiabatic solid inlet face)
    ts_in = np.array([Ts[(ix == 0) & (iy == j)].mean() for j in range(NY)])
    ts_in = (ts_in - TW) / (T0 - TW)
    ts0 = Vs @ c
    print('   entrance |theta_s(col0) - eigen(chi_s=0)| max = %.4f'
          % np.abs(ts_in - np.interp(eta_d, eta, ts0)).max())
    ax = axes[row, 0]
    errs_f, errs_s, errs_x = [], [], []
    for xp, cc in zip(xphys, pcols):
        xs = kappa * xp / (U * A_H**2)
        col = int(round((xp - LX / NX / 2) / (LX / NX)))
        tf = np.array([Tf[(ix == col) & (iy == j)].mean() for j in range(NY)])
        tsg = np.array([Ts[(ix == col) & (iy == j)].mean() for j in range(NY)])
        thf = (tf - TW) / (T0 - TW)
        ths = (tsg - TW) / (T0 - TW)
        af, asol = theta_at(xs, ref)
        af_d = np.interp(eta_d, eta, af)
        as_d = np.interp(eta_d, eta, asol)
        errs_f.append(np.abs(thf - af_d).max())
        errs_s.append(np.abs(ths - as_d).max())
        ax.plot(thf, eta_d, 'o', ms=3, color=cc)
        ax.plot(ths, eta_d, 's', ms=3, mfc='none', color=cc)
        ax.plot(af, eta, '-', color=cc, lw=1, label='x=%gmm' % (xp * 1e3))
        ax.plot(asol, eta, '--', color=cc, lw=1)
        if xp in (0.02, 0.05):     # x-refined (360x40) entrance check
            colx = int(round((xp - LX / 720) / (LX / 360)))
            tfx = np.array([(Tfx[(ixx == colx) & (iyx == j)].mean() - TW)
                            / (T0 - TW) for j in range(NY)])
            ax.plot(tfx, eta_d, '^', ms=3.5, color=cc, mfc='none')
            errs_x.append(np.abs(tfx - af_d).max())
    ax.set_title(title + '  (o fluid, s solid, ^ x-refined; -- solid ref)')
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$\eta=y/a$')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    print('   max profile err at x=20..350 mm: theta_f %.4f theta_s %.4f'
          % (max(errs_f), max(errs_s)))
    print('   x-refined (360x40) entrance theta_f max err %.4f' % max(errs_x))

# fluid bulk temperature vs x* (each run collapses on its own eigenmodes)
axb = axes[0, 1]
for fn, al, col, lab in runs:
    ix, iy, Tf, Ts = load(fn)
    kappa, Pe, _, _, eta, lam, c, Vf, _ = profiles(al)
    cols = np.arange(3, NX - 2, 4)
    xs = kappa * (cols + 0.5) * LX / NX / (U * A_H**2)
    thb = np.array([(Tf[ix == cc].mean() - TW) / (T0 - TW) for cc in cols])
    bulk_mode = np.mean(Vf, axis=0)           # uniform u -> column mean
    xx = np.linspace(0, xs.max() * 1.02, 200)
    thb_a = np.array([np.sum(c * bulk_mode * np.exp(-lam * Pe * v))
                      for v in xx])
    axb.plot(xs, thb, 'o', ms=4, color=col)
    axb.plot(xx, thb_a, '-', color=col, lw=1, label=lab)
axb = axes[1, 1]
for fn, al, col, lab in runs:
    ix, iy, Tf, Ts = load(fn)
    kappa, Pe, _, _, eta, lam, c, Vf, _ = profiles(al)
    cols = np.arange(3, NX - 2, 4)
    xi = (cols + 0.5) * LX / NX / A_H
    thb = np.array([(Tf[ix == cc].mean() - TW) / (T0 - TW) for cc in cols])
    bulk_mode = np.mean(Vf, axis=0)
    xx = np.linspace(0, xi.max() * 1.02, 200)
    thb_a = np.array([np.sum(c * bulk_mode * np.exp(-lam * v))
                      for v in xx])
    axb.plot(xi, thb, 'o', ms=4, color=col)
    axb.plot(xx, thb_a, '-', color=col, lw=1, label=lab)
axb.set_xlabel(r'$\xi=x/a$')
axb.set_ylabel(r'bulk $\theta_{f,b}$')
axb.set_title('Bulk vs physical coordinate')
axb.legend(fontsize=7)
axb.grid(alpha=0.3)
axes[0, 1].set_xlabel(r'$x^*=\kappa_f x/(U a^2)$')
axes[0, 1].set_ylabel(r'bulk $\theta_{f,b}$')
axes[0, 1].set_title('Fluid bulk temperature vs x*')
axes[0, 1].legend(fontsize=7)
axes[0, 1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig('images/ltne_graetz.png', dpi=130)
print('wrote images/ltne_graetz.png')
