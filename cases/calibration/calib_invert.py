"""Inversion (Levenberg-Marquardt) and identifiability (FIM) machinery.

Noise model: each obs group g with relative std s_g contributes
sigma_i = s_g * max(|y_i|, floor_g). Weights w_i = 1/sigma_i^2.

All gradients are computed in the calibration space theta (log10-scaled for
'log' parameters), so the FIM and confidence intervals refer to theta;
relative parameter uncertainty is ln(10)*sigma_theta for log parameters.
"""

import numpy as np

from calib_problems import unpack_theta, theta_true


def group_sigmas(prob, y_ref):
    """Per-observation absolute noise std from the problem's obs_noise list."""
    n = len(y_ref)
    sig = np.empty(n)
    n_grp = len(prob['obs_noise'])
    per = n // n_grp
    for g, (_nm, rel) in enumerate(prob['obs_noise']):
        sl = slice(g * per, (g + 1) * per if g < n_grp - 1 else n)
        floor = max(np.abs(y_ref[sl]).max() * 1e-3, 1e-30)
        sig[sl] = rel * np.maximum(np.abs(y_ref[sl]), floor)
    return sig


def synth_obs(prob, seed=0, theta=None):
    """Synthetic experiment: analytic truth + Gaussian noise."""
    rng = np.random.default_rng(seed)
    th = theta_true(prob) if theta is None else theta
    y_true = prob['forward'](unpack_theta(prob, th))
    sig = group_sigmas(prob, y_true)
    return y_true + rng.normal(0.0, sig), y_true, sig


def jacobian(prob, theta, rel_step=1e-4):
    """Central-difference Jacobian dy/dtheta in calibration space."""
    y0 = prob['forward'](unpack_theta(prob, theta))
    J = np.empty((len(y0), len(theta)))
    for j in range(len(theta)):
        h = rel_step * max(abs(theta[j]), 1.0)
        tp = theta.copy(); tp[j] += h
        tm = theta.copy(); tm[j] -= h
        J[:, j] = (prob['forward'](unpack_theta(prob, tp))
                   - prob['forward'](unpack_theta(prob, tm))) / (2 * h)
    return y0, J


def lm_invert(prob, y_obs, sig, theta0=None, max_it=60, tol=1e-12, verbose=False):
    """Levenberg-Marquardt with box bounds (in theta space)."""
    lo = np.array([np.log10(l) if s == 'log' else l
                   for _n, _t, l, _h, s in prob['params']])
    hi = np.array([np.log10(h) if s == 'log' else h
                   for _n, _t, _l, h, s in prob['params']])
    th = theta0.copy() if theta0 is not None else 0.5 * (lo + hi)
    th = np.clip(th, lo, hi)
    W = 1.0 / sig**2
    lam = 1e-3
    y, J = jacobian(prob, th)
    r = y_obs - y
    cost = np.sum(W * r**2)
    hist = [cost]
    for it in range(max_it):
        JW = J * np.sqrt(W)[:, None]
        A = JW.T @ JW
        g = JW.T @ (r * np.sqrt(W))
        for _ in range(20):
            try:
                d = np.linalg.solve(A + lam * np.diag(np.diag(A)), g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            tn = np.clip(th + d, lo, hi)
            yn = prob['forward'](unpack_theta(prob, tn))
            rn = y_obs - yn
            cn = np.sum(W * rn**2)
            if cn < cost:
                th, y, r, cost = tn, yn, rn, cn
                lam = max(lam / 10, 1e-12)
                break
            lam *= 10
        else:
            break
        hist.append(cost)
        y, J = jacobian(prob, th)
        r = y_obs - y
        if verbose:
            print(f'  it={it:3d} cost={cost:.6e} lam={lam:.1e}')
        if len(hist) > 2 and abs(hist[-2] - hist[-1]) < tol * max(hist[-2], 1e-300):
            break
    return th, np.array(hist)


def fim_analysis(prob, theta):
    """Fisher information analysis at theta.

    Returns dict with normalised sensitivity S_n, FIM, correlation matrix,
    per-parameter relative std (posterior Cramer-Rao bound), condition number.
    """
    y, J = jacobian(prob, theta)
    sig = group_sigmas(prob, y)
    W = 1.0 / sig**2
    FIM = (J * W[:, None]).T @ J
    # normalised sensitivity: d ln y / d ln p (log params) or dy/dp * p/y
    S_n = np.empty_like(J)
    pvals = unpack_theta(prob, theta)
    for j, (nm, _t, _l, _h, s) in enumerate(prob['params']):
        S_n[:, j] = J[:, j] * (np.log(10) if s == 'log' else 1.0 / max(pvals[nm], 1e-300))
    yscale = np.maximum(np.abs(y), np.abs(y).max() * 1e-3)
    S_n = S_n * np.abs(y)[:, None] / yscale[:, None]
    cov = np.linalg.inv(FIM + 1e-300 * np.eye(len(theta)))
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)
    # relative uncertainty of physical params
    rel = np.empty(len(theta))
    for j, (nm, _t, _l, _h, s) in enumerate(prob['params']):
        rel[j] = np.log(10) * sd[j] if s == 'log' else sd[j] / max(abs(pvals[nm]), 1e-300)
    eig = np.linalg.eigvalsh(FIM)
    cond = np.sqrt(eig.max() / max(eig.min(), 1e-300))
    return dict(S_n=S_n, FIM=FIM, corr=corr, rel_std=rel, cond=cond, sig=sig)
