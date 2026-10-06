#!/usr/bin/env python3
# Beavers-Joseph validation case: open 2D channel, upper half clear fluid,
# lower half Darcy porous bed.  Driven by a uniform body force f (N/m^3)
# instead of a pressure gradient, so BOTH ends are pressure-outlet p=0 and
# the driving force per unit volume is exactly f (no dp/dx measurement).
#
#   x in [0,L], L = 0.10 m ; y in [0, 2H], H = 5 mm
#   zone 2 (upper, y in [H,2H]) : fluid
#   zone 3 (lower, y in [0,H])  : porous, VC:porous tag
#
# Analytic (eps=1 so mu_eff=mu, pure Darcy sink mu*u/K in the bed):
#   u_D = K*G/mu,  u_B = (G*H^2/(2*mu) + alpha*sigma*u_D) / (1 + alpha*sigma)
#   sigma = H/sqrt(K), G = f = 0.05 Pa/m
#   fluid profile (y measured from the interface): u = u_B + A*y - G*y^2/(2*mu),
#   A = alpha*(u_B - u_D)/sqrt(K)
#
# Face zones: 4 interior, 5 z-symmetry, 6 bottom wall, 7 top wall,
#             8 outlet p=0 (x=0), 9 outlet p=0 (x=L).
import sys

NX, NY = 100, 40          # NY/2 cells per layer (porous below, fluid above)
LX, LY, LZ = 0.10, 0.010, 0.002
out = sys.argv[1]
MX, MY = NX + 1, NY + 1
NP = MX * MY
NN = 2 * NP
NYC = NY // 2             # cells per zone layer

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
NC = NX * NY              # cells 1..NX*NYC porous (zone 3), rest fluid (zone 2)

L = []
L.append('(0 "generated 2D Beavers-Joseph channel (fluid over porous bed)")')
L.append('(0 "NX=%d NY=%d, zone3=VC:porous lower half, zone2=fluid upper")'
         % (NX, NY))
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
# cell zone 2: fluid (upper half, cell rows j = NYC..NY-1)
NC2 = NX * (NY - NYC)
f2lo, f2hi = NX * NYC + 1, NC
L.append('(0 "Zone 2 %d cells %d..%d, fluid upper half")' % (NC2, f2lo, f2hi))
L.append('(12 (2 %s %s 1 4))' % (hx(f2lo), hx(f2hi)))
L.append('(45 (2 fluid upper)())')
L.append('')
# cell zone 3: porous (lower half, cell rows j = 0..NYC-1)
L.append('(0 "Zone 3 %d cells 1..%d, VC:porous lower half")' % (NX * NYC, NX * NYC))
L.append('(12 (3 1 %s 1 4))' % hx(NX * NYC))
L.append('(45 (3 VC:porous bed)())')
L.append('')

def face(ns, c0, c1):
    L.append(' '.join(hx(v) for v in ns) + ' ' + hx(c0) + ' ' + hx(c1))

# ---- interior faces: y-normal first, then x-normal ----
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
# ---- z symmetry ----
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
# ---- bottom wall y=0 ----
f0, f1 = NF_INT + NF_ZSYM + 1, NF_INT + NF_ZSYM + NF_BOT
L.append('(0 "Zone 6 %d faces %d..%d, wall bottom")' % (NF_BOT, f0, f1))
L.append('(13 (6 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(NX):
    face([bn(i, 0), tn(i, 0), tn(i + 1, 0), bn(i + 1, 0)], cell(i, 0), 0)
L.append('))')
L.append('')
# ---- top wall y=2H ----
f0, f1 = NF_INT + NF_ZSYM + NF_BOT + 1, NF_INT + NF_ZSYM + NF_BOT + NF_TOP
L.append('(0 "Zone 7 %d faces %d..%d, wall top")' % (NF_TOP, f0, f1))
L.append('(13 (7 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(NX):
    face([bn(i, NY), tn(i, NY), tn(i + 1, NY), bn(i + 1, NY)],
         cell(i, NY - 1), 0)
L.append('))')
L.append('')
# ---- outlet x=0 (p=0) ----
f0, f1 = NF_INT + NF_ZSYM + NF_BOT + NF_TOP + 1, \
         NF_INT + NF_ZSYM + NF_BOT + NF_TOP + NF_IN
L.append('(0 "Zone 8 %d faces %d..%d, pressure-outlet x=0")' % (NF_IN, f0, f1))
L.append('(13 (8 %s %s 4 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(0, j), tn(0, j), tn(0, j + 1), bn(0, j + 1)], cell(0, j), 0)
L.append('))')
L.append('')
# ---- outlet x=L (p=0) ----
f0, f1 = NF - NF_OUT + 1, NF
L.append('(0 "Zone 9 %d faces %d..%d, pressure-outlet x=L")' % (NF_OUT, f0, f1))
L.append('(13 (9 %s %s 4 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(NX, j), tn(NX, j), tn(NX, j + 1), bn(NX, j + 1)],
         cell(NX - 1, j), 0)
L.append('))')
L.append('')

with open(out, 'w') as f:
    f.write('\n'.join(L))
print('wrote', out, ': nodes', NN, 'cells', NC, 'faces', NF)
