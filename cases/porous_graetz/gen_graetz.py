#!/usr/bin/env python3
# Graetz forced-convection entry flow in a Darcy porous-filled parallel-plate
# channel.  Only the HALF channel is modelled (symmetry at the centreline):
#   y in [0,a], a = H/2 = 5 mm;  x in [0,L], L = 0.45 m
# Uniform Darcy ("plug") flow along x, inlet T0, constant wall temperature Tw
# at y=a, symmetry (adiabatic) at y=0.  Analytical target is the plug-flow
# Graetz series with effective transverse diffusivity kappa_eff:
#   theta(x,eta) = sum_n 2(-1)^(n+1)/lam_n cos(lam_n eta) exp(-lam_n^2 x*)
#   lam_n = (n-1/2)*pi, eta = y/a, x* = kappa_eff*x/(U*a^2)
#
# Face zones: 3 interior, 4 z-symmetry, 5 centreline symmetry (y=0),
#             6 velocity-inlet (x=0), 7 pressure-outlet (x=L),
#             8 fixed-T wall (y=a).
import sys

NX, NY = 180, 40
LX, LY, LZ = 0.45, 0.005, 0.002
out = sys.argv[1]
MX, MY = NX + 1, NY + 1
NP = MX * MY
NN = 2 * NP

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
NF_CSYM = NX              # centreline y=0
NF_WALL = NX              # heated wall y=a
NF_IN = NY
NF_OUT = NY
NF = NF_INT + NF_ZSYM + NF_CSYM + NF_WALL + NF_IN + NF_OUT
NC = NX * NY

L = []
L.append('(0 "generated 2D porous Graetz channel (half height)")')
L.append('(0 "NX=%d NY=%d, cell zone tagged VC:porous, wall Tw, plug flow")' % (NX, NY))
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
        L.append('  %.15e   %.15e   0.000000000000000e+00' % (i * LX / NX, j * LY / NY))
for j in range(MY):
    for i in range(MX):
        L.append('  %.15e   %.15e   %.15e' % (i * LX / NX, j * LY / NY, LZ))
L.append('))')
L.append('')
L.append('(0 "Zone 2 %d cells 1..%d, VC:porous graetz")' % (NC, NC))
L.append('(12 (2 1 %s 1 4))' % hx(NC))
L.append('(45 (2 VC:porous graetz)())')
L.append('')

def face(ns, c0, c1):
    L.append(' '.join(hx(v) for v in ns) + ' ' + hx(c0) + ' ' + hx(c1))

# ---- interior faces: y-normal first, then x-normal ----
f0, f1 = 1, NF_INT
L.append('(0 "Zone 3 %d faces %d..%d, Interior")' % (NF_INT, f0, f1))
L.append('(13 (3 %s %s 2 4)(' % (hx(f0), hx(f1)))
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
# ---- z symmetry ----
f0, f1 = NF_INT + 1, NF_INT + NF_ZSYM
L.append('(0 "Zone 4 %d faces %d..%d, symmetry")' % (NF_ZSYM, f0, f1))
L.append('(13 (4 %s %s 7 4)(' % (hx(f0), hx(f1)))
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
# ---- centreline symmetry y=0 ----
f0, f1 = NF_INT + NF_ZSYM + 1, NF_INT + NF_ZSYM + NF_CSYM
L.append('(0 "Zone 5 %d faces %d..%d, symmetry centreline")' % (NF_CSYM, f0, f1))
L.append('(13 (5 %s %s 7 4)(' % (hx(f0), hx(f1)))
for i in range(NX):
    face([bn(i, 0), tn(i, 0), tn(i + 1, 0), bn(i + 1, 0)], cell(i, 0), 0)
L.append('))')
L.append('')
# ---- inlet x=0 ----
f0, f1 = NF_INT + NF_ZSYM + NF_CSYM + 1, NF_INT + NF_ZSYM + NF_CSYM + NF_IN
L.append('(0 "Zone 6 %d faces %d..%d, velocity-inlet")' % (NF_IN, f0, f1))
L.append('(13 (6 %s %s 5 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(0, j), tn(0, j), tn(0, j + 1), bn(0, j + 1)], cell(0, j), 0)
L.append('))')
L.append('')
# ---- outlet x=L ----
f0, f1 = NF_INT + NF_ZSYM + NF_CSYM + NF_IN + 1, NF_INT + NF_ZSYM + NF_CSYM + NF_IN + NF_OUT
L.append('(0 "Zone 7 %d faces %d..%d, pressure-outlet")' % (NF_OUT, f0, f1))
L.append('(13 (7 %s %s 4 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(NX, j), tn(NX, j), tn(NX, j + 1), bn(NX, j + 1)],
         cell(NX - 1, j), 0)
L.append('))')
L.append('')
# ---- heated wall y=a ----
f0, f1 = NF - NF_WALL + 1, NF
L.append('(0 "Zone 8 %d faces %d..%d, wall fixed-T")' % (NF_WALL, f0, f1))
L.append('(13 (8 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(NX):
    face([bn(i, NY), tn(i, NY), tn(i + 1, NY), bn(i + 1, NY)],
         cell(i, NY - 1), 0)
L.append('))')
L.append('')

with open(out, 'w') as f:
    f.write('\n'.join(L))
print('wrote', out, ': nodes', NN, 'cells', NC, 'faces', NF)
