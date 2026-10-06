#!/usr/bin/env python3
"""Nusselt-number validation charts: computed values (parsed from solver
logs) vs literature references, for the pure-fluid and porous benchmarks."""
import re, sys, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def nu_from_log(fn):
    m = re.findall(r'Nusselt number.*?:\s*([0-9.eE+-]+)', open(fn).read())
    return float(m[-1]) if m else None

def main(logdir, out):
    # ---- pure-fluid Rayleigh-Benard (hot bottom) --------------------------
    fluid = {}   # Ra -> Nu  (40x40 unless noted)
    for ra, lg in [(1e3,'rb_Ra1000.log'), (3e3,'rb_Ra3000.log'),
                   (5e3,'rb_Ra5000.log'), (1e4,'rb_Ra1e4.log'),
                   (5e4,'rb_Ra5e4.log')]:
        fluid[ra] = nu_from_log(f'{logdir}/fluid/{lg}')
    fluid_80 = nu_from_log(f'{logdir}/fluid/rb_Ra1e4_80.log')
    # Hollands (1976) infinite horizontal air layer, Pr=0.71 anchors
    hollands = {5e3: 1.90, 1e4: 2.39, 5e4: 3.44}

    # ---- porous side-heated Darcy cavity ----------------------------------
    por40 = {10: nu_from_log(f'{logdir}/porous_darcy/porous_Ra10.log'),
             100: nu_from_log(f'{logdir}/porous_darcy/porous_Ra100.log'),
             1000: nu_from_log(f'{logdir}/porous_darcy/verify_Ra1000.log')}
    por80 = nu_from_log(f'{logdir}/porous_darcy/por_Ra1000_80.log')
    por160 = nu_from_log(f'{logdir}/porous_darcy/verify_Ra1000_160.log')
    lit = {10: [(1.07,'Walker & Homsy')],
           100: [(3.10,'Bejan'), (3.16,'Baytas & Pop')],
           1000: [(13.64,'Mahmud & Fraser'), (14.06,'Baytas & Pop')]}

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6))

    a = ax[0]
    ras = sorted(fluid)
    a.plot(ras, [fluid[r] for r in ras], 'o-', color='tab:red', lw=1.5,
           ms=7, label='computed (40x40)')
    a.plot([1e4], [fluid_80], 's', color='tab:red', mfc='w', ms=8,
           label='computed (80x80)')
    hr = sorted(hollands)
    a.plot(hr, [hollands[r] for r in hr], 'D--', color='k', mfc='none',
           lw=1.2, ms=7, label='Hollands (1976), infinite layer')
    a.axvline(1708, color='gray', ls=':', lw=1)
    a.text(1750, 1.05, 'onset $Ra=1708$', rotation=90, fontsize=8,
           color='gray', va='bottom')
    a.set_xscale('log'); a.set_xlabel('Rayleigh number  Ra')
    a.set_ylabel('Nusselt number  Nu'); a.set_ylim(0.5, 4.2)
    a.set_title('Pure-fluid Rayleigh-Benard (heated from below)')
    a.legend(fontsize=8, loc='upper left'); a.grid(alpha=0.3)

    b = ax[1]
    pr = sorted(por40)
    b.plot(pr, [por40[r] for r in pr], 'o-', color='tab:blue', lw=1.5, ms=7,
           label='computed (40x40)')
    b.plot([1000, 1000], [por40[1000], por160], '-', color='tab:blue', lw=0.8)
    b.plot([1000], [por80], 'd', color='tab:blue', ms=6, label='computed (80x80)')
    b.plot([1000], [por160], 's', color='tab:blue', mfc='w', ms=8,
           label='computed (160x160)')
    for ra, pts in lit.items():
        for nu, who in pts:
            b.plot([ra], [nu], '^', color='k', mfc='none', ms=8)
            b.annotate(who, (ra, nu), textcoords='offset points', xytext=(6, -3),
                       fontsize=7)
    b.plot([], [], '^', color='k', mfc='none', ms=8, label='literature')
    b.axvline(4*np.pi**2, color='gray', ls=':', lw=1)
    b.text(4*np.pi**2*1.05, 1.05, 'onset $4\\pi^2$', rotation=90,
           fontsize=8, color='gray', va='bottom')
    b.set_xscale('log'); b.set_xlabel('Darcy-Rayleigh number  $Ra_K$')
    b.set_ylabel('Nusselt number  Nu'); b.set_ylim(0.5, 16)
    b.set_title('Porous Darcy cavity (side-heated, LTE)')
    b.legend(fontsize=8, loc='upper left'); b.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print('wrote', out)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
