"""Calibration problem definitions for porous-media parameter identification.

Each problem is a dict with keys:
  name        str
  params      list of (name, true_value, lower, upper, scale)
              scale is 'log' or 'lin' (log -> calibrate in log10 space)
  obs_points  measurement configuration (problem-specific)
  forward(p)  reduced-order analytic model: p dict -> obs vector (1-D ndarray)
  obs_noise   relative std of synthetic measurement noise per obs group
  case_ref    existing case directory providing solver forward runs

Reduced-order models duplicate the semi-analytic references already validated
in cases/porous_plug, cases/porous_disp, cases/ltne_disp, cases/beavers_joseph.
The solver VTU provides "truth" data (with discretisation error); the analytic
model is the inversion engine (fast). This avoids the inverse crime of using
the same discrete model for both synthesis and inversion.
"""

import numpy as np
from math import erf as _erf_scalar

_erf = np.vectorize(_erf_scalar)

# --------------------------------------------------------------------------
# Shared fluid / porous constants (300 K air, matching validated cases)
# --------------------------------------------------------------------------
RHO = 1.177          # kg/m3
MU = 1.846e-5        # Pa.s
CP = 1005.0          # J/kg/K
KF = 0.026           # W/m/K (fluid conductivity)


# ==========================================================================
# P1: permeability K + Forchheimer coefficient C_F  (porous-plug 1-D)
# ==========================================================================
# Pressure drop of a 1-D plug, length L, superficial velocity u:
#   dp = L * ( mu/(eps*K) * u + rho*C_F*u^2 )
# Measure dp at several velocities -> quadratic in u separates K and C_F.

P1_L = 0.1           # m  (plug length, cases/porous_plug)
P1_EPS = 0.4

def _p1_forward(p):
    K = p['K']; CF = p['C_F']
    u = P1_PROB['obs_points']['u']
    return P1_L * (MU / (P1_EPS * K) * u + RHO * CF * u**2)

P1_PROB = dict(
    name='P1_plug_K_CF',
    params=[('K',   1.0e-8, 1e-10, 1e-6, 'log'),
            ('C_F', 500.0,  1.0,   5e3,  'log')],
    obs_points=dict(
        # 6 velocities spanning Darcy -> Forchheimer dominated regimes
        # Re_p = rho u sqrt(K)/mu : 0.032 .. 12.8
        u=np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0]),
    ),
    forward=_p1_forward,
    obs_noise=[('dp', 0.01)],          # 1% pressure transducer
    case_ref='../porous_plug',
)


# ==========================================================================
# P2: transverse dispersion alpha_t  (porous-disp mixing layer)
# ==========================================================================
# Temperature profile at axial station x0:
#   T(y) = T_mid - dT/2 * erf( (y-H/2) * sqrt(u / (4*kappa*x0)) )
#   kappa = (k_cond + rho cp alpha_t u) / (rho cp)
# Half-width of the erf grows with sqrt(alpha_t).

P2_U = 1.0
P2_H = 0.06
P2_TMID = 305.0
P2_DT = 50.0          # large hot/cold step: TC noise scales with |T| (~0.6 K),
                      # so the profile signal must be >> 1 K for identifiability

def _p2_forward(p):
    at = p['alpha_t']
    cfg = P2_PROB['obs_points']
    out = []
    for x0 in cfg['x_stations']:
        kappa = (KF + RHO * CP * at * P2_U) / (RHO * CP)
        eta = (cfg['y'] - P2_H / 2) * np.sqrt(P2_U / (4 * kappa * x0))
        out.append(P2_TMID - P2_DT / 2 * _erf(eta))
    return np.concatenate(out)

P2_PROB = dict(
    name='P2_mixing_alpha_t',
    params=[('alpha_t', 1.0e-3, 1e-5, 1e-1, 'log')],
    obs_points=dict(
        # three axial stations, 15 transverse points each
        x_stations=[0.01, 0.025, 0.04],
        y=np.linspace(0.015, 0.045, 15),
    ),
    forward=_p2_forward,
    obs_noise=[('T', 0.002)],          # 0.2% (~0.6 K on 305 K) thermocouple
    case_ref='../porous_disp',
)


# ==========================================================================
# P3: longitudinal dispersion alpha_l + interstitial h_sf (LTNE sweat cooling)
# ==========================================================================
# 1-D coupled ODE (identical to cases/ltne_disp/plot_ltne_disp.py reference):
#   G cp Tf' = (Kf Tf')' + H (Ts - Tf)      Kf = eps*kf + rho cp alpha_l u
#   0        = Ds Ts''   + H (Tf - Ts)      Ds = (1-eps) ks
# BCs: Tf(0)=Tin, Ds Ts'(0)=0, Kf Tf'(L)=qf, Ds Ts'(L)=qs
# q split by microscopic conduction fraction (Nield).

P3_L = 0.05
P3_EPS = 0.3
P3_KS = 16.2
P3_G = 2.0           # kg/m2/s
P3_TIN = 300.0
P3_Q = 2.0e6         # W/m2

def _p3_solve(alpha_l, H, xsamp):
    """Semi-analytic coupled-ODE solution sampled at xsamp (from inlet)."""
    eps = P3_EPS; L = P3_L
    a = P3_G * CP
    u = P3_G / RHO
    Kf = eps * KF + RHO * CP * alpha_l * u
    Ds = (1 - eps) * P3_KS
    # cubic characteristic roots  [Ds*Kf, -Ds*a, -H(Ds+Kf), H a]
    coeff = [Ds * Kf, -Ds * a, -H * (Ds + Kf), H * a]
    r = np.roots(coeff)
    r = r[np.argsort(r.real)]
    # unknowns [c0, d1, d2, d3] with d_j = c_j exp(r_j L) (stable shift)
    # Tf(x) = c0 + sum d_j exp(r_j (x-L));  Ts = m_j * modal amplitude
    m = 1.0 + (a * r - Kf * r**2) / H
    q = P3_Q
    qf = q * Kf / (Kf + Ds)
    qs = q * Ds / (Kf + Ds)
    A = np.zeros((4, 4), complex)
    b = np.zeros(4, complex)
    # Tf(0) = Tin
    A[0, 0] = 1.0
    A[0, 1:] = np.exp(-r * L)
    b[0] = P3_TIN
    # Ds Ts'(0) = 0
    A[1, 1:] = m * r * np.exp(-r * L)
    # Kf Tf'(L) = qf
    A[2, 1:] = Kf * r
    b[2] = qf
    # Ds Ts'(L) = qs
    A[3, 1:] = Ds * m * r
    b[3] = qs
    sol = np.linalg.solve(A, b)
    c0 = sol[0].real
    d = sol[1:]
    E = np.exp(np.outer(xsamp - L, r))
    Tf = c0 + E @ d
    Ts = E @ (m * d)
    return Tf.real, Ts.real

def _p3_forward(p):
    al = p['alpha_l']; H = p['h_sf']
    xs = P3_PROB['obs_points']['x']
    Tf, Ts = _p3_solve(al, H, xs)
    return np.concatenate([Tf, Ts])

P3_PROB = dict(
    name='P3_ltne_alpha_l_hsf',
    params=[('alpha_l', 5.0e-3, 1e-4, 5e-2, 'log'),
            ('h_sf',    2.0e5,  1e3,   1e7,  'log')],
    obs_points=dict(
        # 3 interior + outlet stations; both phases measured
        x=np.array([0.015, 0.035, 0.05]),
    ),
    forward=_p3_forward,
    obs_noise=[('Tf', 0.005), ('Ts', 0.005)],   # 0.5% high-temp TC
    case_ref='../ltne_disp',
)


# ==========================================================================
# P4: Beavers-Joseph slip coefficient bj_alpha (open channel over porous bed)
# ==========================================================================
# Two-layer stress-continuity + velocity-jump model (4x4 coupled ODE), the
# continuum limit of the discrete BJ flux validated in cases/beavers_joseph.
# Channel: fluid half 0<y<Hf (eps=1), porous half -Hp<y<0 with perm K.
# Body force f drives both. Interface at y=0 with slip coefficient alpha.

P4_HF = 0.01         # fluid half-height
P4_HP = 0.01         # porous half-height
P4_K = 1.0e-6
P4_F = 0.05          # body force N/m3
P4_EPS = 1.0         # validated case used eps=1 (Brinkman mu_eff=mu)

def _p4_forward(p):
    alpha = p['bj_alpha']
    Ny = P4_PROB['obs_points']['ny']
    # solve the 4x4 coupled ODE on a fine grid by finite differences, then
    # sample.  Model (non-dimensional y in [-Hp, Hf]):
    #   fluid : mu u_f'' + f = 0
    #   porous: mu u_p'' - mu/K u_p + f = 0
    #   interface: u_f - u_p = (sqrt(K)/alpha) (u_f' + u_p')/2  (BJ Robin)
    #              u_f'(0) = u_p'(0)                            (stress cont.)
    #   walls: u_f(Hf)=0, u_p(-Hp)=0
    # Direct analytic solution is available; use it.
    mu = MU; K = P4_K; f = P4_F
    lam = np.sqrt(K)
    Hf = P4_HF; Hp = P4_HP
    # porous particular u_pp = f K / mu; general u_p = A cosh(y/lam) + B sinh(y/lam) + fK/mu
    # fluid u_f = -f/(2mu) y^2 + C y + D
    # BCs:
    #  u_f(Hf) = 0
    #  u_p(-Hp) = 0
    #  stress: u_f'(0) = u_p'(0)
    #  BJ: u_f(0) - u_p(0) = (lam/alpha) * u_f'(0)     (slip, single-sided grad)
    # (The discrete validated model uses stress continuity + velocity jump;
    #  see cases/beavers_joseph/README.md. We use the same two conditions.)
    upp = f * K / mu
    # unknowns A, B, C, D
    M = np.zeros((4, 4)); rhs = np.zeros(4)
    # u_f(Hf)=0
    M[0, 2] = Hf; M[0, 3] = 1.0
    rhs[0] = f * Hf**2 / (2 * mu)
    # u_p(-Hp)=0
    M[1, 0] = np.cosh(-Hp / lam); M[1, 1] = np.sinh(-Hp / lam)
    rhs[1] = -upp
    # stress continuity at 0: C = B/lam
    M[2, 1] = -1.0 / lam; M[2, 2] = 1.0
    # BJ jump: D - (A + upp) = (lam/alpha) C
    M[3, 0] = -1.0; M[3, 2] = -lam / alpha; M[3, 3] = 1.0
    rhs[3] = upp
    A, B, C, D = np.linalg.solve(M, rhs)
    y = np.concatenate([np.linspace(-Hp, 0, Ny // 2 + 1)[:-1],
                        np.linspace(0, Hf, Ny - Ny // 2)])
    u = np.where(y <= 0,
                 A * np.cosh(y / lam) + B * np.sinh(y / lam) + upp,
                 -f / (2 * mu) * y**2 + C * y + D)
    return u

P4_PROB = dict(
    name='P4_bj_alpha',
    params=[('bj_alpha', 1.0, 0.05, 10.0, 'log')],
    obs_points=dict(ny=21),      # PIV-like vertical velocity profile
    forward=_p4_forward,
    obs_noise=[('u', 0.02)],     # 2% PIV
    case_ref='../beavers_joseph',
)


PROBLEMS = {p['name']: p for p in (P1_PROB, P2_PROB, P3_PROB, P4_PROB)}


def pack_theta(prob, values):
    """Physical parameter values -> calibration vector theta (log-scaled)."""
    th = []
    for (name, _tv, _lo, _hi, scale), v in zip(prob['params'], values):
        th.append(np.log10(v) if scale == 'log' else v)
    return np.array(th)


def unpack_theta(prob, theta):
    """Calibration vector -> dict of physical parameter values."""
    p = {}
    for (name, _tv, lo, hi, scale), t in zip(prob['params'], theta):
        v = 10**t if scale == 'log' else t
        p[name] = np.clip(v, lo, hi)
    return p


def theta_true(prob):
    return pack_theta(prob, [tv for _n, tv, _l, _h, _s in prob['params']])
