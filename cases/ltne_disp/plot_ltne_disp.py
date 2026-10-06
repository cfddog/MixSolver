#!/usr/bin/env python3
# LTNE 1D slab + longitudinal thermal dispersion: numerical T_f / T_s vs
# the EXACT steady two-phase ODE solution (4th order, three exponential
# modes + constant, coefficients from a 4x4 BC linear system).
#
# Per-unit-bulk-volume steady equations (u = superficial velocity = G/rho):
#   G*cp*T_f' = Kf*T_f'' + H*(T_s - T_f)
#   0         = Ds*T_s'' + H*(T_f - T_s)
# with
#   Kf = eps*k_f + rho*cp*disp_l*u    (fluid axial conduction + dispersion)
#   Ds = (1-eps)*k_s                  (solid conduction)
#   H  = h_sf*a_sf
#
# Eliminate T_s -> 4th-order ODE for T_f:
#   Ds*Kf*Tf'''' - Ds*a*Tf''' - H*(Ds+Kf)*Tf'' + H*a*Tf' = 0,  a = G*cp
# modes: Tf = A0 + sum_j Aj*exp(rj*x), rj = 3 roots of the cubic;
# Ts follows each mode with multiplier ms(r) = 1 + (a*r - Kf*r^2)/H.
#
# BCs matching the FVM implementation exactly:
#   Tf(0) = Tin                    (Dirichlet fluid inlet)
#   Ds*Ts'(0) = 0                  (solid adiabatic at inlet faces)
#   Kf*Tf'(L) = qf, Ds*Ts'(L) = qs (outlet flux split, microscopic
#                                   conductivities only, dispersion
#                                   excluded from the split as in code)
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- case parameters (must match the .control files) ----------------------
RHO, MU, CP, KF0 = 1.177, 1.846e-5, 1005.0, 0.026
EPS, KS = 0.3, 16.2
G, Q, L, TIN = 2.0, 2.0e6, 0.05, 300.0
HA = 2.0e5
NX = 100
U = G / RHO
DS = (1 - EPS) * KS
A_CONV = G * CP
KFT = EPS * KF0                       # microscopic fluid conductivity
# ---------------------------------------------------------------------------


def reference(disp_l):
    Kf = KFT + RHO * CP * disp_l * U
    qf = Q * KFT / (KFT + DS)          # code split ignores dispersion
    qs = Q * DS / (KFT + DS)
    roots = np.roots([DS * Kf, -DS * A_CONV, -HA * (DS + Kf), HA * A_CONV])
    roots = np.real_if_close(roots, tol=1000)
    if np.iscomplexobj(roots) and np.abs(roots.imag).max() > 1e-8:
        raise RuntimeError('unexpected complex roots: %s' % roots)
    r = np.real(roots)
    ms = 1.0 + (A_CONV * r - Kf * r**2) / HA

    # Unknowns [c0, d1, d2, d3] with dj = cj*exp(rj*L): keeps every exponent
    # bounded on x in [0,L] (stiff fluid-conduction mode, esp. disp_l=0).
    em = np.exp(-r * L)
    M = np.zeros((4, 4))
    b = np.zeros(4)
    # (1) Tf(0) = Tin:  c0 + sum dj*exp(-rj*L) = Tin
    M[0, 0] = 1.0
    M[0, 1:] = em
    b[0] = TIN
    # (2) Ds*Ts'(0) = 0
    M[1, 1:] = ms * r * em
    # (3) Kf*Tf'(L) = qf
    M[2, 1:] = Kf * r
    b[2] = qf
    # (4) Ds*Ts'(L) = qs
    M[3, 1:] = DS * ms * r
    b[3] = qs
    coef = np.linalg.solve(M, b)
    c0, d = coef[0], coef[1:]

    def Tf(x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        return c0 + np.sum(d[:, None] * np.exp(r[:, None] * (x[None, :] - L)),
                           axis=0)

    def Ts(x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        return c0 + np.sum(d[:, None] * ms[:, None]
                           * np.exp(r[:, None] * (x[None, :] - L)), axis=0)

    return Tf, Ts, Kf, qf, qs


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
    Tf = arr('temperature', cd_)
    Ts = arr('temperature_solid', cd_)
    pts = arr('Coordinates', re.search(r'<Points>(.*?)</Points>', s, re.S)
              .group(1)).reshape(-1, 3)
    st = np.concatenate(([0], off[:-1]))
    dx = L / NX
    ix = np.zeros(len(off), dtype=int)
    for i in range(len(off)):               # hex cells exported as 12 tets:
        ids = conn[st[i]:off[i]]            # strip leading nverts token
        if len(ids) == 5 and ids[0] == 4:
            ids = ids[1:]
        ix[i] = int(round((pts[ids, 0].mean() - dx / 2) / dx))
    xl = (np.arange(NX) + 0.5) * dx
    Tfl = np.array([Tf[ix == i].mean() for i in range(NX)])
    Tsl = np.array([Ts[ix == i].mean() for i in range(NX)])
    return xl, Tfl, Tsl


runs = [('ltne_disp0.vtu', 0.0, 'C0', 'baseline (disp_l=0)'),
        ('ltne_disp1.vtu', 0.005, 'C3', 'disp_l=0.005 m')]

print('%12s %9s %10s %10s %10s %10s %10s'
      % ('run', 'Kf W/mK', 'dTf(L)n', 'dTf(L)a', 'errTf mx', 'errTf md',
         'errTs mx'))
data = {}
for fn, al, col, lab in runs:
    xl, Tfn, Tsn = load(fn)
    Tfa, Tsa, Kf, qf, qs = reference(al)
    data[al] = (xl, Tfn, Tsn, Tfa, Tsa)
    rise = max(abs(Tfa(np.array([L]))[0] - TIN), 1.0)
    ef = np.abs(Tfn - Tfa(xl)) / rise * 100
    es = np.abs(Tsn - Tsa(xl)) / (Tsa(xl) - TIN).max() * 100
    print('%12s %9.3f %10.2f %10.2f %9.2f%% %9.2f%% %9.2f%%'
          % (lab.split()[0], Kf, Tfn[-1] - TIN,
             Tfa(np.array([L]))[0] - TIN, ef.max(), np.median(ef), es.max()))

xa = np.linspace(0, L, 400)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for fn, al, col, lab in runs:
    xl, Tfn, Tsn, Tfa, Tsa = data[al]
    ax[0].plot(xa, Tfa(xa), '-', color=col, lw=1.3)
    ax[0].plot(xl, Tfn, 'o', ms=3.5, mfc='none', color=col, label=lab)
    ax[1].plot(xa, Tsa(xa), '-', color=col, lw=1.3)
    ax[1].plot(xl, Tsn, 'o', ms=3.5, mfc='none', color=col, label=lab)
    ax[2].plot(xa, Tsa(xa) - Tfa(xa), '-', color=col, lw=1.3)
    ax[2].plot(xl, Tsn - Tfn, 'o', ms=3.5, mfc='none', color=col, label=lab)

ax[0].set_xlabel('x (m)'); ax[0].set_ylabel('$T_f$ [K]')
ax[0].set_title('Fluid temperature')
ax[1].set_xlabel('x (m)'); ax[1].set_ylabel('$T_s$ [K]')
ax[1].set_title('Solid temperature')
ax[2].set_xlabel('x (m)'); ax[2].set_ylabel(r'$T_s-T_f$ [K]')
ax[2].set_title('Phase temperature difference')
for a in ax:
    a.grid(alpha=0.3); a.legend(fontsize=8)
fig.suptitle(r'LTNE 1D slab: longitudinal dispersion (Bear) enters the '
             r'fluid phase only;  $G=2$ kg/m$^2$s, $q''=2e6$ W/m$^2$')
fig.tight_layout()
fig.savefig('images/ltne_disp.png', dpi=130)
print('wrote images/ltne_disp.png')
