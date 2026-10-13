#!/usr/bin/env python3
"""Cylinder Re = 200 (UNSTEADY) scene check: is the wake shedding, and at what
Strouhal number?

Usage:
    python3 check_cylinder_unsteady.py <dt> <log|snapshot> [<...> ...]

  <dt>        the PIMPLE time step used in unMesh.control (e.g. 0.005).  It is
              the only quantity the script cannot recover from the files.
  <log>       logs/<tag>.log -- the per-step solver log.  PIMPLE prints one line
              per time step ("step=  N  t= ...  du/dt= ...  max|u|= ..."), so
              this gives the FINEST possible time series.
  <snapshot>  out/<tag>.step<NNNN>.vtu -- one transient VTU snapshot per
              n_out_every steps.  Coarser but device-independent.

You can pass a log, snapshots, or both; each source is analysed on its own.

Background
  At Re = 200 the 2D wake is a von Karman street (periodic limit cycle,
  St = f*D/U ~ 0.20, mean CD ~ 1.35).  A steady solver returns the symmetric
  base flow instead, so the tests are
    * "are the signals periodic?"  -> the antisymmetric (shedding) mode is on
    * "what is the period?"        -> St = f*D/U, expect ~0.2 (T ~ 5 at U=D=1)

What is printed for each source
  * span, number of samples, sample interval, detrended RMS and peak-to-peak
  * the dominant frequency (windowed FFT over the developed part of the record)
    and the resulting Strouhal number
  * for a VTU snapshot series: u_y at a wake probe (x = 3, y = +/-0.6) and the
    whole-field max|u_y|, which also grows when the wake rolls up

If the signals come out flat (rms ~ 0) the run is the steady symmetric solution
and the shedding has not started yet -- extend n_time_max, or see the shedding
notes in README.md.
"""
import re
import sys

import numpy as np

D_CYL = 1.0          # cylinder diameter
U_INF = 1.0          # free-stream speed
X_PROBE = 3.0        # wake probe: 2.5 D behind the rear stagnation point
Y_PROBE = 0.6        # off-axis, where the transverse fluctuation is largest


# --------------------------------------------------------------------------
# signal analysis
# --------------------------------------------------------------------------
def analyse(t, y, label, unit=''):
    """Print detrended amplitude + dominant frequency (Strouhal if unit-less)."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 8:
        print('  %-22s : only %d samples -- skipping' % (label, n))
        return
    dt = float(np.median(np.diff(t)))
    # Statistics on the END of the record (the start is the growth transient):
    # the last 40%, but never fewer than 16 samples.
    n_use = max(16, int(0.4 * n))
    i0 = max(0, n - n_use)
    ts, ys = t[i0:], y[i0:]
    # linear detrend, then remove the mean of the detrended signal
    coef = np.polyfit(ts, ys, 1)
    yl = ys - np.polyval(coef, ts)
    yl = yl - yl.mean()
    rms = float(np.sqrt(np.mean(yl ** 2)))
    ptp = float(yl.max() - yl.min())
    print('  %-22s : samples %d, dt %.4g, span %.4g'
          % (label, n, dt, t[-1] - t[0]))
    print('    mean %.5e   detrended rms %.3e   peak-to-peak %.3e %s'
          % (float(ys.mean()), rms, ptp, unit))

    # windowed FFT on the detrended signal
    m = len(yl)
    w = np.hanning(m)
    Y = np.fft.rfft(yl * w)
    fr = np.fft.rfftfreq(m, dt)
    if len(fr) < 3:
        return
    k = int(np.argmax(np.abs(Y[1:])) + 1)
    # parabolic interpolation on |Y| around the peak -> sub-bin frequency.
    # This matters: the raw bin spacing 1/(window span) would otherwise be the
    # dominant uncertainty on the Strouhal number.
    dk = 0.0
    if 0 < k < len(Y) - 1:
        y0, y1, y2 = np.abs(Y[k - 1]), np.abs(Y[k]), np.abs(Y[k + 1])
        den = y0 - 2.0 * y1 + y2
        if abs(den) > 0.0:
            dk = 0.5 * (y0 - y2) / den
            dk = max(-0.5, min(0.5, dk))
    df = float(fr[1] - fr[0]) if len(fr) > 1 else 0.0
    f = float(fr[k]) + dk * df
    amp = 2.0 * float(np.abs(Y[k])) / max(w.sum(), 1e-30)
    print('    window %.4g (df %.4g), peak bin %d @ %.5g, interp shift %+.2f'
          % (ts[-1] - ts[0], df, k, fr[k], dk))
    print('    dominant frequency %.5g 1/time  (period %.5g, amplitude %.3e)'
          % (f, 1.0 / f if f > 0 else np.inf, amp))
    if unit == '':
        print('    -> Strouhal number St = f*D/U = %.4f   (reference ~0.20 '
              'at Re=200)' % (f * D_CYL / U_INF))
    if rms < 1e-8:
        print('    -> FLAT: steady symmetric solution, no shedding yet')
    else:
        print('    -> oscillatory: the antisymmetric shedding mode is active')


def read_log(fn, dt):
    """Parse 'step=  N  t= ...  du/dt= ...  max|u|= ...' lines (one per step)."""
    t, dudt, umax = [], [], []
    pat = re.compile(r'step=\s*(\d+)\s+t=\s*([-\d.eE+]+).*?'
                     r'du/dt=\s*([-\d.eE+]+).*?max\|u\|=\s*([-\d.eE+]+)')
    for line in open(fn, errors='replace'):
        if 'Tmin=' in line:
            continue
        m = pat.search(line)
        if m:
            t.append(float(m.group(2)))
            dudt.append(float(m.group(3)))
            umax.append(float(m.group(4)))
    print('log %s : %d step records' % (fn, len(t)))
    if not t:
        print('  no "step=" records found -- is this a transient log?')
        return
    analyse(t, umax, 'max|u| (log)')
    analyse(t, dudt, 'du/dt (log)')


# --------------------------------------------------------------------------
# VTU snapshot series
# --------------------------------------------------------------------------
def read_vtu_vel(fn):
    """Cell-centroid coordinates + velocity from a bin/mixnsolver .vtu."""
    txt = open(fn).read()

    def arr(name, cols=1, dtype=float):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % re.escape(name),
                      txt, re.S)
        if m is None:
            raise SystemExit('check_cylinder_unsteady: no "%s" array in %s'
                             % (name, fn))
        a = np.array([dtype(x) for x in m.group(1).split()])
        return a.reshape(-1, cols) if cols > 1 else a

    pts = arr('Coordinates', 3)
    con = arr('connectivity', 5, int)
    vel = arr('velocity', 3)
    cen_id = con[:, 4]
    first = np.concatenate(([True], cen_id[1:] != cen_id[:-1]))
    idx = np.nonzero(first)[0]
    return pts[cen_id[idx], :], vel[idx, :]


def probe_velocity(xyz, vel, xq, knn=8):
    """IDW over the knn nearest cell centroids (as mod_uns_output probe_velocity)."""
    d2 = np.sum((xyz - xq) ** 2, axis=1)
    idx = np.argsort(d2, kind='stable')[:knn]
    w = 1.0 / (d2[idx] + 1.0e-12)
    return (w[:, None] * vel[idx]).sum(axis=0) / w.sum()


def read_snapshots(files, dt):
    files = sorted(files,
                   key=lambda f: int(re.search(r'step(\d+)', f).group(1)))
    t, vy_p, vy_m, umax, vymax = [], [], [], [], []
    zp = None
    for fn in files:
        step = int(re.search(r'step(\d+)', fn).group(1))
        xyz, vel = read_vtu_vel(fn)
        if zp is None:
            zp = 0.5 * (xyz[:, 2].min() + xyz[:, 2].max())
            print('snapshots: %d files, mesh cells %d, z-probe %.3f'
                  % (len(files), len(xyz), zp))
        t.append(step * dt)
        vy_p.append(float(probe_velocity(xyz, vel,
                                         np.array([X_PROBE, +Y_PROBE, zp]))[1]))
        vy_m.append(float(probe_velocity(xyz, vel,
                                         np.array([X_PROBE, -Y_PROBE, zp]))[1]))
        umax.append(float(np.hypot(vel[:, 0], vel[:, 1]).max()))
        vymax.append(float(np.abs(vel[:, 1]).max()))
    if len(t) > 1:
        print('  (snapshot interval %.4g, first t = %.3f, last t = %.3f)'
              % (t[1] - t[0], t[0], t[-1]))
    analyse(t, np.array(vy_p), 'u_y(x=3, y=+0.6)')
    analyse(t, np.array(umax), 'max|u| (field)')
    analyse(t, np.array(vy_p) - np.array(vy_m), 'u_y(+0.6)-u_y(-0.6)')
    if t:
        print('  last snapshot: max|u_y| = %.5f, |u_y(+0.6)-u_y(-0.6)| = %.5f'
              % (vymax[-1], abs(vy_p[-1] - vy_m[-1])))


# --------------------------------------------------------------------------
def main(argv):
    if len(argv) < 3:
        raise SystemExit('usage: check_cylinder_unsteady.py <dt> '
                         '<log|snapshot.vtu> [...]')
    dt = float(argv[1])
    # .vtu (and .vtu_stepNNNNNNN.vtu) are snapshots, anything else is a log
    vtus = [a for a in argv[2:] if '.vtu' in a]
    logs = [a for a in argv[2:] if '.vtu' not in a]
    if not logs and not vtus:
        raise SystemExit('check_cylinder_unsteady: no log or snapshot given')
    for fn in logs:
        print('')
        read_log(fn, dt)
    if vtus:
        print('')
        read_snapshots(vtus, dt)


if __name__ == '__main__':
    main(sys.argv)
