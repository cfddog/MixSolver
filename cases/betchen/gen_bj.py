#!/usr/bin/env python3
# Betchen 2006 validation case 1 (BJ): fluid channel over a porous bed,
# interface PARALLEL to the flow, uniform-velocity inlet + pressure outlet.
#
#   x in [0, L], L = 8H ; y in [0, 2H], H the half channel height.
#   Lower half y in [0,H]  : porous (zone 3, VC:porous)
#   Upper half y in [H,2H] : fluid   (zone 2)
#
# Interface condition (13)/(14): velocity + (tangential) stress continuity,
# which the current solver reproduces with the DEFAULT internal face
# (no bj_alpha).  Driven by a uniform inlet velocity U0 (Beta: Re_H=1):
#   Re_H = rho*U0*H/mu = 1  =>  U0 = mu/(rho*H)
# Outlet: pressure-outlet p=0.  Walls y=0 and y=2H no-slip; z symmetry.
#
# Physical parameters (paper Sec 5.1):
#   eps = 0.7, Da = K/H^2 = 1e-2 or 1e-3, cE = 1.75*eps/(150*eps^5)^0.5
import sys

H = 0.01                 # half channel height (m)
LX, LY, LZ = 8 * H, 2 * H, 0.002   # L=8H, 2H full height, 2 z-layers
NX, NY = 100, 40         # paper grid ny=40, nx=100
out = sys.argv[1]
MX, MY = NX + 1, NY + 1
NP = MX * MY
NN = 2 * NP
NYC = NY // 2            # cells per zone layer (porous lower, fluid upper)
rho, mu = 1.177, 1.846e-5
U0 = mu / (rho * H)      # Re_H = 1


def bn(i, j):
    return 1 + i + j * MX


def tn(i, j):
    return bn(i, j) + NP


def cell(i, j):
    return 1 + i + j * NX


def hx(v):
    return format(v, 'x')


NF_YINT = (NY - 1) * NX
NF_XINT = (NX - 1) * NY
NF_INT = NF_YINT + NF_XINT
NF_ZSYM = 2 * NX * NY
NF_BOT = NX
NF_TOP = NX
NF_IN = NY
NF_OUT = NY
NF = NF_INT + NF_ZSYM + NF_BOT + NF_TOP + NF_IN + NF_OUT
NC = NX * NY

L = []
L.append('(0 "generated Betchen BJ channel (fluid over porous bed)")')
L.append('(0 "NX=%d NY=%d, zone3=VC:porous lower half, zone2=fluid upper")'
         % (NX, NY))
L.append('(0 "H=%.4f m L=8H=%.4f m, Re_H=1 -> U0=%.6e m/s")' % (H, LX, U0))
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
    for i in range(MX):
        L.append('  %.15e   %.15e   0.000000000000000e+00'
                 % (i * LX / NX, j * LY / NY))
for j in range(MY):
    for i in range(MX):
        L.append('  %.15e   %.15e   %.15e' % (i * LX / NX, j * LY / NY, LZ))
L.append('))')
L.append('')
NC2 = NX * (NY - NYC)
f2lo, f2hi = NX * NYC + 1, NC
L.append('(0 "Zone 2 %d cells %d..%d, fluid upper half")' % (NC2, f2lo, f2hi))
L.append('(12 (2 %s %s 1 4))' % (hx(f2lo), hx(f2hi)))
L.append('(45 (2 fluid upper)())')
L.append('')
L.append('(0 "Zone 3 %d cells 1..%d, VC:porous lower half")'
         % (NX * NYC, NX * NYC))
L.append('(12 (3 1 %s 1 4))' % hx(NX * NYC))
L.append('(45 (3 VC:porous bed)())')
L.append('')


def face(ns, c0, c1):
    L.append(' '.join(hx(v) for v in ns) + ' ' + hx(c0) + ' ' + hx(c1))


f0, f1 = 1, NF_INT
L.append('(0 "Zone 4 %d faces %d..%d, Interior")' % (NF_INT, f0, f1))
L.append('(13 (4 %s %s 2 4)(' % (hx(f0), hx(f1)))
for j in range(1, NY):
    for i in range(NX):
        face([bn(i, j), tn(i, j), tn(i + 1, j), bn(i + 1, j)],
             cell(i, j - 1), cell(i, j))
for i in range(1, NX):
    for j in range(NY):
        face([bn(i, j), tn(i, j), tn(i, j + 1), bn(i, j + 1)],
             cell(i - 1, j), cell(i, j))
L.append('))')
L.append('')
f0, f1 = NF_INT + 1, NF_INT + NF_ZSYM
L.append('(0 "Zone 5 %d faces %d..%d, symmetry")' % (NF_ZSYM, f0, f1))
L.append('(13 (5 %s %s 7 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    for i in range(NX):
        face([bn(i, j), bn(i + 1, j), bn(i + 1, j + 1), bn(i, j + 1)],
             cell(i, j), 0)
for j in range(NY):
    for i in range(NX):
        face([tn(i, j), tn(i + 1, j), tn(i + 1, j + 1), tn(i, j + 1)],
             cell(i, j), 0)
L.append('))')
L.append('')
f0, f1 = NF_INT + NF_ZSYM + 1, NF_INT + NF_ZSYM + NF_BOT
L.append('(0 "Zone 6 %d faces %d..%d, wall bottom")' % (NF_BOT, f0, f1))
L.append('(13 (6 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(NX):
    face([bn(i, 0), tn(i, 0), tn(i + 1, 0), bn(i + 1, 0)], cell(i, 0), 0)
L.append('))')
L.append('')
f0, f1 = NF_INT + NF_ZSYM + NF_BOT + 1, NF_INT + NF_ZSYM + NF_BOT + NF_TOP
L.append('(0 "Zone 7 %d faces %d..%d, wall top")' % (NF_TOP, f0, f1))
L.append('(13 (7 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(NX):
    face([bn(i, NY), tn(i, NY), tn(i + 1, NY), bn(i + 1, NY)],
         cell(i, NY - 1), 0)
L.append('))')
L.append('')
f0, f1 = NF_INT + NF_ZSYM + NF_BOT + NF_TOP + 1, \
         NF_INT + NF_ZSYM + NF_BOT + NF_TOP + NF_IN
L.append('(0 "Zone 8 %d faces %d..%d, velocity-inlet x=0")' % (NF_IN, f0, f1))
L.append('(13 (8 %s %s 5 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(0, j), tn(0, j), tn(0, j + 1), bn(0, j + 1)], cell(0, j), 0)
L.append('))')
L.append('')
f0, f1 = NF - NF_OUT + 1, NF
L.append('(0 "Zone 9 %d faces %d..%d, pressure-outlet x=L")'
         % (NF_OUT, f0, f1))
L.append('(13 (9 %s %s 4 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(NX, j), tn(NX, j), tn(NX, j + 1), bn(NX, j + 1)],
         cell(NX - 1, j), 0)
L.append('))')
L.append('')

with open(out, 'w') as f:
    f.write('\n'.join(L))
print('wrote', out, ': nodes', NN, 'cells', NC,
      'faces', NF, 'U0=%.6e m/s' % U0)