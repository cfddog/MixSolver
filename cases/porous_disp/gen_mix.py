#!/usr/bin/env python3
# 2D porous thermal mixing layer: uniform Darcy flow along x through a
# VC:porous block; the inlet is split at mid-height into a hot (T1) lower half
# and a cold (T2) upper half.  Transverse thermal dispersion (alpha_t) widens
# the mixing layer, giving the error-function similarity solution
#   T(x,y) = (T1+T2)/2 + (T2-T1)/2 * erf( (y-H/2) * sqrt(u/(4*Gamma_y*x)) )
#   Gamma_y = k_cond + rho*cp*alpha_t*u     (u = superficial velocity)
#
# Face zones: 3 interior, 4 z-symmetry, 5 y-symmetry,
#             6 velocity-inlet lower half (T1), 7 pressure-outlet,
#             8 velocity-inlet upper half (T2).
import sys

NX, NY = 80, 80
LX, LY, LZ = 0.04, 0.06, 0.002
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
NF_YSYM = 2 * NX
NF_IN1 = NY // 2          # lower half inlet
NF_IN2 = NY // 2          # upper half inlet
NF_OUT = NY
NF = NF_INT + NF_ZSYM + NF_YSYM + NF_IN1 + NF_IN2 + NF_OUT
NC = NX * NY

L = []
L.append('(0 "generated 2D porous thermal mixing layer")')
L.append('(0 "NX=%d NY=%d, cell zone tagged VC:porous, split T inlet")' % (NX, NY))
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
L.append('(0 "Zone 2 %d cells 1..%d, VC:porous mixing layer")' % (NC, NC))
L.append('(12 (2 1 %s 1 4))' % hx(NC))
L.append('(45 (2 VC:porous mixlayer)())')
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
# ---- y symmetry (top and bottom walls of the channel) ----
f0, f1 = NF_INT + NF_ZSYM + 1, NF_INT + NF_ZSYM + NF_YSYM
L.append('(0 "Zone 5 %d faces %d..%d, symmetry")' % (NF_YSYM, f0, f1))
L.append('(13 (5 %s %s 7 4)(' % (hx(f0), hx(f1)))
for i in range(NX):   # y=0
    face([bn(i, 0), tn(i, 0), tn(i + 1, 0), bn(i + 1, 0)], cell(i, 0), 0)
for i in range(NX):   # y=H
    face([bn(i, NY), tn(i, NY), tn(i + 1, NY), bn(i + 1, NY)], cell(i, NY - 1), 0)
L.append('))')
L.append('')
# ---- inlet lower half (hot) ----
f0, f1 = NF_INT + NF_ZSYM + NF_YSYM + 1, NF_INT + NF_ZSYM + NF_YSYM + NF_IN1
L.append('(0 "Zone 6 %d faces %d..%d, velocity-inlet T1")' % (NF_IN1, f0, f1))
L.append('(13 (6 %s %s 5 4)(' % (hx(f0), hx(f1)))
for j in range(NY // 2):
    face([bn(0, j), tn(0, j), tn(0, j + 1), bn(0, j + 1)], cell(0, j), 0)
L.append('))')
L.append('')
# ---- outlet ----
f0, f1 = NF_INT + NF_ZSYM + NF_YSYM + NF_IN1 + NF_IN2 + 1, NF
L.append('(0 "Zone 7 %d faces %d..%d, pressure-outlet")' % (NF_OUT, f0, f1))
L.append('(13 (7 %s %s 4 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(NX, j), tn(NX, j), tn(NX, j + 1), bn(NX, j + 1)],
         cell(NX - 1, j), 0)
L.append('))')
L.append('')
# ---- inlet upper half (cold) ----
f0, f1 = NF_INT + NF_ZSYM + NF_YSYM + NF_IN1 + 1, \
         NF_INT + NF_ZSYM + NF_YSYM + NF_IN1 + NF_IN2
L.append('(0 "Zone 8 %d faces %d..%d, velocity-inlet T2")' % (NF_IN2, f0, f1))
L.append('(13 (8 %s %s 5 4)(' % (hx(f0), hx(f1)))
for j in range(NY // 2, NY):
    face([bn(0, j), tn(0, j), tn(0, j + 1), bn(0, j + 1)], cell(0, j), 0)
L.append('))')
L.append('')

with open(out, 'w') as f:
    f.write('\n'.join(L))
print('wrote', out, ': nodes', NN, 'cells', NC, 'faces', NF)
