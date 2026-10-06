"""Driver: synthetic calibration closure for all problems.

Usage: python3 run_calib.py [P1 P2 P3 P4]   (default: all)

For each problem:
  1. synthetic obs = analytic truth + noise (seeded, reproducible)
  2. LM inversion from mid-bounds start -> recovered params + rel. error
  3. FIM analysis at the truth -> correlation matrix, Cramer-Rao bounds
  4. figure: forward fit vs obs (+ FIM correlation inset info in stdout)
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from calib_problems import PROBLEMS, unpack_theta, theta_true
from calib_invert import synth_obs, lm_invert, fim_analysis


def run(name, ax_fit):
    prob = PROBLEMS[name]
    pnames = [p[0] for p in prob['params']]
    y_obs, y_true, sig = synth_obs(prob, seed=1)

    # deterministic off-truth start: +0.3 dex (factor ~2) on every log param,
    # clipped to bounds -- avoids trivial starts at the true value
    th0 = theta_true(prob) + 0.3
    lo = np.array([np.log10(l) if s == 'log' else l
                   for _n, _t, l, _h, s in prob['params']])
    hi = np.array([np.log10(h) if s == 'log' else h
                   for _n, _t, _l, h, s in prob['params']])
    th0 = np.clip(th0, lo, hi)
    th_hat, hist = lm_invert(prob, y_obs, sig, theta0=th0)
    p_hat = unpack_theta(prob, th_hat)
    p_true = unpack_theta(prob, theta_true(prob))

    fa = fim_analysis(prob, theta_true(prob))

    print(f'\n=== {name} ===')
    print(f'  LM converged: cost {hist[0]:.3e} -> {hist[-1]:.3e} in {len(hist)-1} it')
    hdr = f'  {"param":10s} {"true":>12s} {"recovered":>12s} {"rel.err":>9s} {"CRB rel.std":>11s}'
    print(hdr)
    for j, nm in enumerate(pnames):
        tv = p_true[nm]; hv = p_hat[nm]
        print(f'  {nm:10s} {tv:12.4e} {hv:12.4e} {abs(hv-tv)/tv:9.2e} {fa["rel_std"][j]:11.2e}')
    print(f'  FIM cond = {fa["cond"]:.3e}')
    if len(pnames) > 1:
        print('  correlation matrix:')
        for j, nm in enumerate(pnames):
            row = ' '.join(f'{fa["corr"][j,k]:+7.3f}' for k in range(len(pnames)))
            print(f'    {nm:10s} {row}')

    y_fit = prob['forward'](p_hat)
    n_grp = len(prob['obs_noise'])
    per = len(y_obs) // n_grp
    for g in range(n_grp):
        sl = slice(g * per, (g + 1) * per if g < n_grp - 1 else len(y_obs))
        x = np.arange(sl.start, sl.stop)
        ax_fit.errorbar(x, y_obs[sl], yerr=sig[sl], fmt='o', ms=3, lw=0.8,
                        label=f'obs grp{g}')
        ax_fit.plot(x, y_true[sl], '-', lw=1.2, label=f'truth grp{g}')
        ax_fit.plot(x, y_fit[sl], '--', lw=1.0, label=f'fit grp{g}')
    ax_fit.set_title(name)
    ax_fit.set_xlabel('obs index'); ax_fit.set_ylabel('measurement')
    ax_fit.legend(fontsize=7)
    return hist


def main():
    sel = sys.argv[1:] or list(PROBLEMS)
    names = [n for n in PROBLEMS if any(s in n for s in sel)]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, name in zip(axes.ravel(), names):
        run(name, ax)
    fig.tight_layout()
    fig.savefig('images/calib_closure.png', dpi=140)
    print('\nfigure: images/calib_closure.png')


if __name__ == '__main__':
    main()
