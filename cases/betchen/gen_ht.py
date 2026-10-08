#!/usr/bin/env python3
# Betchen 2006 validation case 3 (HT): bottom-heated aluminium-foam block,
# LTNE two-temperature.  Case A: foam block only.  Case B: block plus a 2 mm
# upstream pure-air gap (fluid/porous conjugate, tests the energy matching).
#
#   block : y in [0,H], x in [xg, xg+L]   (H=45 mm, L=114 mm)
#   gap   : x in [0, xg]          (Case B only, xg=2 mm)
#   inlet at x=0 (velocity-inlet), outlet x=xg+L (pressure-outlet)
#   bottom y=0 heated (tbc), top y=H adiabatic wall, z symmetry
#
# Foam (Calmidi-Mahajan Sample 2): eps=0.9118, K=1.8e-7, cE=0.085,
#   k_se=6.46, k_fe=0.0237.  Air at 300 K.  h_sf from Nu_sf correlation.
#
# Cell zones:
#   Case A : zone 2 porous (whole block)
#   Case B : zone 2 fluid (gap) + zone 3 porous (block)
# Cells must be x-ordered so each zone is one contiguous id range (x slowest).
# Face zones: 5 interior, 6 z-sym, 7 wall top (y=H), 8 wall bottom (y=0),
#             9 velocity-inlet, a pressure-outlet.
import sys

H_ht = 0.045                            # block height (y, mm->m)
LB = 0.114                              # block length (x, mm->m)
EPS = 0.9118
K = 1.8e-7
RHO, MU = 1.177, 1.846e-5
CP, KF = 1005.0, 0.026
U = 1.0                                 # mean through velocity
N_Y = 70                                # control volumes in y (paper grid)
N_X_BLOCK = 80                          # in x over the block
N_X_GAP = 7                             # gap CVs (Case B), paper configure
A = "block"
if len(sys.argv) > 2 and sys.argv[2] == 'gap':
    A = "gap"
                                    # xg = 2 mm only for Case B
XG = 0.002 if A == "gap" else 0.0
L_TOT = LB + XG
N_X = N_X_BLOCK + (N_X_GAP if A == "gap" else 0)
out = sys.argv[1]

LZ = 0.002
MX, MY = N_X + 1, N_Y + 1
NP = MX * MY
NN = 2 * NP

# --- Calmidi-Mahajan h_sf correlation --------------------------------
def nu_sf(re_dl, pr):
    CT = 0.52
    return CT * re_dl ** 0.5 * pr ** 0.37

d_l = 0.55e-3
pr = MU * CP / KF
re_dl = RHO * U * d_l / MU
h_sf = nu_sf(re_dl, pr) * KF / d_l
# volumetric surface area of the metal foam (approx, m^2/m^3)
d_p = 3.8e-3
a_sf = (1 - EPS) * 6.0 / d_p
HVOL = h_sf * a_sf

def bn(i, j):
    return 1 + i + j * MX

def tn(i, j):
    return bn(i, j) + NP


def xnode(i):
    # Case A: uniform over [0,L]. Case B: gap (0..XG) split into N_X_GAP,
    # then block (XG..XG+L) split into N_X_BLOCK.
    if A == "block":
        return LB * i / N_X
    if i <= N_X_GAP:
        return XG * i / N_X_GAP
    return XG + LB * (i - N_X_GAP) / N_X_BLOCK


def cell(i, j):
    # x slowest so gap zone (i<N_X_GAP) and block zone are contiguous ranges
    return 1 + j + i * N_Y


def hx(v):
    return format(v, 'x')


NF_YINT = (N_Y - 1) * N_X
NF_XINT = (N_X - 1) * N_Y
NF_INT = NF_YINT + NF_XINT
NF_ZSYM = 2 * N_X * N_Y
NF_WTOP = N_X
NF_WBOT = N_X
NF_IN = N_Y
NF_OUT = N_Y
NF = NF_INT + NF_ZSYM + NF_WTOP + NF_WBOT + NF_IN + NF_OUT
NC = N_X * N_Y

L = []
L.append('(0 "generated Betchen HT aluminium-foam block (LTNE)")')
L.append('(0 "N_X=%d N_Y=%d case=%s, L=%.3f H=%.3f")'
         % (N_X, N_Y, A, LB, H_ht))
L.append('(0 "h_sf=%.1f W/m2K a_sf=%.1f /m h*a=%.1f W/m3K Re_dl=%.1f")'
         % (h_sf, a_sf, HVOL, re_dl))
L.append('')
L.append('(0 "Dimension : 3")')
L.append('(2 3)')
L.append('')
L.append('(0 "Number of Nodes : %d")' % NN)
L.append('(10 (0 1 %s 0 3))' % hx(NN))
L.append('')
L.append('(0 "Total Number of Faces : %d")' % NF)
L.append('(13 (0 1 %s 0))' % hx(NF))
L.append('')
L.append('(0 "Total Number of Cells : %d")' % NC)
L.append('(12 (0 1 %s 0))' % hx(NC))
L.append('')
L.append('(0 "Zone 1  Number of Nodes : %d")' % NN)
L.append('(10 (1 1 %s 1 3)(' % hx(NN))
for j in range(MY):
    yj = H_ht * j / N_Y
    for i in range(MX):
        L.append('  %.15e   %.15e   0.000000000000000e+00'
                 % (xnode(i), yj))
for j in range(MY):
    yj = H_ht * j / N_Y
    for i in range(MX):
        L.append('  %.15e   %.15e   %.15e' % (xnode(i), yj, LZ))
L.append('))')
L.append('')
# ---- cell zones ------------------------------------------------------
if A == "block":
    z2_lo = cell(0, 0)
    z2_hi = cell(N_X - 1, N_Y - 1)
    L.append('(0 "Zone 2 %d cells %d..%d, VC:porous foam block")'
             % (NC, z2_lo, z2_hi))
    L.append('(12 (2 %s %s 1 4))' % (hx(z2_lo), hx(z2_hi)))
    L.append('(45 (2 VC:porous block)())')
    L.append('')
else:
    g_hi = cell(N_X_GAP - 1, N_Y - 1)
    b_lo = g_hi + 1
    b_hi = cell(N_X - 1, N_Y - 1)
    L.append('(0 "Zone 2 %d cells %d..%d, fluid air gap")'
             % (N_X_GAP * N_Y, 1, g_hi))
    L.append('(12 (2 %s %s 1 4))' % (hx(1), hx(g_hi)))
    L.append('(45 (2 fluid gap)())')
    L.append('')
    L.append('(0 "Zone 3 %d cells %d..%d, VC:porous foam block")'
             % (N_X_BLOCK * N_Y, b_lo, b_hi))
    L.append('(12 (3 %s %s 1 4))' % (hx(b_lo), hx(b_hi)))
    L.append('(45 (3 VC:porous block)())')
    L.append('')


def face(ns, c0, c1):
    L.append(' '.join(hx(v) for v in ns) + ' ' + hx(c0) + ' ' + hx(c1))


f0, f1 = 1, NF_INT
L.append('(0 "Zone 5 %d faces %d..%d, Interior")' % (NF_INT, f0, f1))
L.append('(13 (5 %s %s 2 4)(' % (hx(f0), hx(f1)))
for j in range(1, N_Y):
    for i in range(N_X):
        face([bn(i, j), tn(i, j), tn(i + 1, j), bn(i + 1, j)],
             cell(i, j - 1), cell(i, j))
for i in range(1, N_X):
    for j in range(N_Y):
        face([bn(i, j), tn(i, j), tn(i, j + 1), bn(i, j + 1)],
             cell(i - 1, j), cell(i, j))
L.append('))')
L.append('')
f0, f1 = NF_INT + 1, NF_INT + NF_ZSYM
L.append('(0 "Zone 6 %d faces %d..%d, symmetry")' % (NF_ZSYM, f0, f1))
L.append('(13 (6 %s %s 7 4)(' % (hx(f0), hx(f1)))
for j in range(N_Y):
    for i in range(N_X):
        face([bn(i, j), bn(i + 1, j), bn(i + 1, j + 1), bn(i, j + 1)],
             cell(i, j), 0)
for j in range(N_Y):
    for i in range(N_X):
        face([tn(i, j), tn(i + 1, j), tn(i + 1, j + 1), tn(i, j + 1)],
             cell(i, j), 0)
L.append('))')
L.append('')
# top wall y=H
f0, f1 = NF_INT + NF_ZSYM + 1, NF_INT + NF_ZSYM + NF_WTOP
L.append('(0 "Zone 7 %d faces %d..%d, wall top y=H")' % (NF_WTOP, f0, f1))
L.append('(13 (7 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(N_X):
    face([bn(i, N_Y), tn(i, N_Y), tn(i + 1, N_Y), bn(i + 1, N_Y)],
         cell(i, N_Y - 1), 0)
L.append('))')
L.append('')
# bottom wall y=0 (heated)
f0, f1 = NF_INT + NF_ZSYM + NF_WTOP + 1, NF_INT + NF_ZSYM + NF_WTOP + NF_WBOT
L.append('(0 "Zone 8 %d faces %d..%d, wall bottom y=0 (heated)")'
         % (NF_WBOT, f0, f1))
L.append('(13 (8 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(N_X):
    face([bn(i, 0), tn(i, 0), tn(i + 1, 0), bn(i + 1, 0)], cell(i, 0), 0)
L.append('))')
L.append('')
# inlet x=0
f0, f1 = NF_INT + NF_ZSYM + NF_WTOP + NF_WBOT + 1, \
         NF_INT + NF_ZSYM + NF_WTOP + NF_WBOT + NF_IN
L.append('(0 "Zone 9 %d faces %d..%d, velocity-inlet x=0")'
         % (NF_IN, f0, f1))
L.append('(13 (9 %s %s 5 4)(' % (hx(f0), hx(f1)))
for j in range(N_Y):
    face([bn(0, j), tn(0, j), tn(0, j + 1), bn(0, j + 1)], cell(0, j), 0)
L.append('))')
L.append('')
# outlet x=L
f0, f1 = NF - NF_OUT + 1, NF
L.append('(0 "Zone a %d faces %d..%d, pressure-outlet x=L")'
         % (NF_OUT, f0, f1))
# zone id 'a' = hex 10 (see gen_plug.py note)
L.append('(13 (a %s %s 4 4)(' % (hx(f0), hx(f1)))
for j in range(N_Y):
    face([bn(N_X, j), tn(N_X, j), tn(N_X, j + 1), bn(N_X, j + 1)],
         cell(N_X - 1, j), 0)
L.append('))')
L.append('')

with open(out, 'w') as f:
    f.write('\n'.join(L))
print('wrote', out, ': case=%s nodes=%d cells=%d faces=%d'
      % (A, NN, NC, NF))
print('  h_sf=%.1f a_sf=%.1f h*a=%.1f Re_dl=%.1f Pr=%.3f'
      % (h_sf, a_sf, HVOL, re_dl, pr))