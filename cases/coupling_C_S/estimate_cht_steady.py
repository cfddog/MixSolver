#!/usr/bin/env python3
"""
estimate_cht_steady.py -- analytical anchor for the coupling_C_S CHT case.

Computes, from the case's own control files and mesh data, the laminar
flat-plate anchor that the coupled run must approach:

  * the model fluid the OpenCFD-EC nondimensional solver actually sees: it
    fixes mu* = 1/Re plus gamma/Pr, hence cp = gamma*R/(gamma-1) and
    k = mu*cp/Pr (NOT the real-air k of the reference state);
  * free-stream state and Blasius boundary layer thickness;
  * the adiabatic wall (recovery) temperature T_aw;
  * the coupled steady state of the fluid-film + slab series pair
        T_w = (h*T_aw + (k_s/L_s)*T_base) / (h + k_s/L_s)
        q   = h*(T_aw - T_w) = (k_s/L_s)*(T_w - T_base)
    in three flavours -- average h, local h at x = L, and the first-cell
    conductance k_f/d_1 -- which bracket the Dirichlet-Neumann loop gain
    G = h*R_s that iface_t_relax must damp;
  * the recommended fixed relaxation w = 1/(1+G) for each bound.

Optionally recover the *measured* loop gain from two consecutive interface
residuals of a run (the '[struct rank 0] diag iter ... CHT relax: w, |r|max'
lines): r_{k+1}/r_k = 1 - w*(1+G).  A *negative* recovered gain means the
residual growth cannot be explained by the relaxation feedback at all, i.e.
the drift is physical -- then changing w is pointless.  (Measured on
coupling_C_S: G_meas = -21.5 from iterations 4 -> 5 with w = 0.005.)

USAGE
-----
    python3 estimate_cht_steady.py [--json]
    python3 estimate_cht_steady.py --measure-resid 102.0 141.0 --relax 0.03
"""

from __future__ import annotations

import argparse
import json
import math
import sys

# ---- case inputs (mirror control.ec / mix.control / unMesh.control) -------
R_GAS = 287.058          # J/kg/K
GAMMA = 1.4
PR_L = 0.7               # control.ec PrL
MA = 3.0                 # control.ec Ma
T_INF = 800.0            # control.ec T_inf             [K]
T_REF = 800.0            # mix.control T_ref            [K]
RHO_REF = 0.44124        # mix.control rho_ref          [kg/m^3]
U_REF = 1700.95          # mix.control u_ref = Ma*a_ref [m/s]
L_REF = 0.001            # mix.control L_ref = 1 mm     [m]
RE_GRID = 1000.0         # control.ec Re (per grid unit, i.e. per L_REF)

L_PLATE = 1.0            # fluid plate length = interface length [m]
A_IFACE = 0.1            # interface area 1000 mm x 100 mm     [m^2]
D_FIRST = 1.0e-6         # struct first cell height (check_iface_alignment) [m]
K_S = 16.2               # unMesh.control cell_zone k_s  [W/m/K]
CP_S = 500.0
RHO_S = 8000.0
L_S = 0.1                # slab thickness (y -100..0 mm) [m]
T_BASE = 300.0           # unMesh.control tbc zone 4     [K]


def main() -> int:
    ap = argparse.ArgumentParser(description="coupling_C_S analytical anchor")
    ap.add_argument("--measure-resid", nargs=2, type=float, metavar=("R0", "R1"),
                    help="two consecutive |r|max interface residuals [K]")
    ap.add_argument("--relax", type=float, default=None,
                    help="the iface_t_relax used for --measure-resid")
    ap.add_argument("--json", action="store_true", help="also dump json")
    args = ap.parse_args()

    ok = True

    # ---------------- reference state ----------------
    a_ref = math.sqrt(GAMMA * R_GAS * T_REF)
    u_from_ma = MA * a_ref
    p_inf = RHO_REF * R_GAS * T_REF
    cp = GAMMA * R_GAS / (GAMMA - 1.0)
    print("--- reference state (mix.control) ---")
    print(f"  a_ref    = sqrt(gamma R T_ref) = {a_ref:.2f} m/s")
    print(f"  Ma*a_ref = {u_from_ma:.2f} m/s   (u_ref = {U_REF:.2f} m/s)")
    if abs(u_from_ma - U_REF) / U_REF > 1e-3:
        print("  WARN u_ref does not match Ma*a_ref (>0.1%): the struct's "
              "nondimensional velocity scale is the free stream, they must agree")
        ok = False
    print(f"  p_inf    = rho_ref R T_ref = {p_inf:.1f} Pa "
          f"(1 atm = 101325; struct p* = 1/(gamma Ma^2))")
    if abs(p_inf - 101325.0) / 101325.0 > 0.01:
        print("  WARN p_inf differs from 1 atm by more than 1%")
        ok = False

    # ---------------- model fluid seen by the struct solver ----------------
    mu = RHO_REF * U_REF * L_REF / RE_GRID
    nu = mu / RHO_REF
    k_f = mu * cp / PR_L
    print("--- model fluid (nondim solver: mu* = 1/Re, cp = gamma R/(gamma-1)) ---")
    print(f"  mu  = rho_ref*u_ref*L_ref/Re = {mu:.6e} Pa.s "
          f"(air at 800 K is 3.7e-5 -> viscosity-enlarged model fluid)")
    print(f"  nu  = {nu:.6e} m^2/s")
    print(f"  cp  = gamma R/(gamma-1)      = {cp:.1f} J/kg/K")
    print(f"  k_f = mu*cp/Pr               = {k_f:.4f} W/m/K")

    # ---------------- boundary layer ----------------
    re_l = RHO_REF * U_REF * L_PLATE / mu
    delta = 5.0 * L_PLATE / math.sqrt(re_l)
    print("--- laminar flat plate (Blasius) ---")
    print(f"  Re_L (L = 1 m) = {re_l:.3e}   Re_x(1 mm) = {re_l * L_REF:.3f}")
    print(f"  delta(x=1 m) = 5.0*x/sqrt(Re_x) = {delta * 1e3:.2f} mm; "
          f"first cell {D_FIRST * 1e6:.1f} um -> {delta / D_FIRST:.0f} cells "
          f"across the layer")
    h_avg = 0.664 * math.sqrt(re_l) * PR_L ** (1.0 / 3.0) * k_f / L_PLATE
    h_loc = 0.332 * math.sqrt(re_l) * PR_L ** (1.0 / 3.0) * k_f / L_PLATE
    h_cell = k_f / D_FIRST
    print(f"  h_avg  (0.664*sqrt(Re)*Pr^1/3*k/L) = {h_avg:.1f} W/m^2/K")
    print(f"  h_local(x = L)                     = {h_loc:.1f} W/m^2/K")
    print(f"  k_f/d_1 (first-cell conductance)   = {h_cell:.4e} W/m^2/K")

    # ---------------- adiabatic wall temperature ----------------
    print("--- adiabatic wall (recovery) temperature ---")
    t_aw = {}
    for name, r in (("Pr^1/3 (laminar, used below)", PR_L ** (1.0 / 3.0)),
                    ("sqrt(Pr)", math.sqrt(PR_L)),
                    ("1.0 (total temperature)", 1.0)):
        t_aw[name] = T_INF * (1.0 + r * 0.5 * (GAMMA - 1.0) * MA ** 2)
        print(f"  r = {name:28s} = {r:.4f} -> T_aw = {t_aw[name]:.1f} K")
    t_aw_use = t_aw["Pr^1/3 (laminar, used below)"]

    # ---------------- coupled steady state ----------------
    r_s = L_S / K_S
    print("--- coupled steady state (series: fluid film + slab) ---")
    print(f"  R_s = L_s/k_s = {r_s:.4e} m^2K/W   "
          f"(alpha_s = {K_S / (RHO_S * CP_S):.3e} m^2/s)")
    rows = []
    for name, h in (("h_avg", h_avg), ("h_local", h_loc), ("k_f/d_1", h_cell)):
        t_w = (h * t_aw_use + (K_S / L_S) * T_BASE) / (h + K_S / L_S)
        q = h * (t_aw_use - t_w)
        g = h * r_s
        rows.append(dict(case=name, h=h, T_w=t_w, q=q, Q=q * A_IFACE,
                         G=g, w_opt=1.0 / (1.0 + g)))
        print(f"  {name:8s} h = {h:10.3e} -> T_w = {t_w:7.1f} K, "
              f"q = {q:10.3e} W/m^2, Q = {q * A_IFACE:9.3e} W, "
              f"G = {g:9.2f}, w_opt = {1.0 / (1.0 + g):.3e}")
    print(f"  (T_aw = {t_aw_use:.1f} K, T_base = {T_BASE:.0f} K, "
          f"A_iface = {A_IFACE} m^2)")
    print("  --> the Dirichlet-Neumann gain G lies between the h_local and the")
    print("      k_f/d_1 row; only a run settles which one applies.  Watch, per")
    print("      coupling iteration, '[uns rank 1] diag ... Twall min/max' (must")
    print("      climb from 300 K towards T_w) and 'CHT energy: q_iface_in,")
    print("      q_base_in, resid' (resid should shrink against q_iface_in).")

    # ---------------- optional measured gain ----------------
    if args.measure_resid:
        r0, r1 = args.measure_resid
        if not args.relax:
            print("WARN --measure-resid needs --relax <w> to recover G")
        else:
            w = args.relax
            ratio = r1 / r0 if r0 else float("nan")
            g_meas = (1.0 - ratio) / w - 1.0
            print("--- measured loop gain ---")
            print(f"  |r| {r0:g} -> {r1:g} K, ratio = {ratio:.4f}, w = {w:g}")
            print(f"  G_meas = (1 - ratio)/w - 1 = {g_meas:.3f}")
            if g_meas > 0.0:
                print(f"  suggested w = 0.7/(1 + G) = {0.7 / (1.0 + g_meas):.4e} "
                      f"(0.7 = 30% stability margin)")
            else:
                print("  G_meas <= 0: the residual growth is NOT explained by the")
                print("  T-relaxation feedback (a loop gain G > -1 would give")
                print("  |r| ratio < 1 + w*(1+G) ~ 1+w).  The drift is therefore")
                print("  the physical transient, not a relaxation instability --")
                print("  changing w will not fix it (see cases/coupling_C_S/README.md")
                print("  section 7.4/8: a cold slab soaking up a ~1e8 W/m^2")
                print("  transient interface flux).")

    if args.json:
        print(json.dumps(dict(a_ref=a_ref, p_inf=p_inf, mu=mu, nu=nu, cp=cp,
                              k_f=k_f, Re_L=re_l, delta=delta, h_avg=h_avg,
                              h_loc=h_loc, h_cell=h_cell, T_aw=t_aw, R_s=r_s,
                              steady=rows), indent=2))

    print("RESULT: " + ("OK" if ok else "WARN (see messages above)"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
