#!/usr/bin/env python3
# Generate a 1D porous-plug hex .cas mesh: NX x NY x 1 wedge-thick duct.
# Coordinates are already in metres.  The single cell zone carries a
# space-free "VC:porous" tag in its (45 record, so the solver must infer the
# porous block type from the zone name (the .control file supplies the porous
# coefficients via a type-less cell_zone line).
#
# Face zones: 3 interior, 4 z-symmetry, 5 y-symmetry,
#             6 velocity-inlet (x=0), 7 pressure-outlet (x=L).
# All Fluent indices are hexadecimal.
import sys

NX, NY = 40, 2
LX, LY, LZ = 0.1, 0.01, 0.002
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
NF_IN = NY
NF_OUT = NY
NF = NF_INT + NF_ZSYM + NF_YSYM + NF_IN + NF_OUT
NC = NX * NY

L = []
L.append('(0 "generated 1D porous plug")')
L.append('(0 "NX=%d NY=%d, cell zone tagged VC:porous")' % (NX, NY))
L.append('')
L.append('(0 "Dimension : 3")')
L.append('(2 3)')
L.append('')
L.append('(0 "Number of Nodes : %d")' % NN)
L.append('(10 (0 1 %s 0 3))' % hx(NN))
L.append('')
L.append('(0 "Total Number of Faces : %d")' % NF)
L.append('(0 "       Boundary Faces : %d")' % (NF_ZSYM + NF_YSYM + NF_IN + NF_OUT))
L.append('(0 "       Interior Faces : %d")' % NF_INT)
L.append('(13 (0 1 %s 0))' % hx(NF))
L.append('')
L.append('(0 "Total Number of Cells : %d")' % NC)
L.append('(0 "            Hex cells : %d")' % NC)
L.append('(12 (0 1 %s 0))' % hx(NC))
L.append('')
# ---- nodes ----
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
# ---- cells (uniform hex, no payload) ----
L.append('(0 "Zone 2 %d cells 1..%d, VC:porous plug")' % (NC, NC))
L.append('(12 (2 1 %s 1 4))' % hx(NC))
# (45 record: cond-name token holds the space-free VC tag (the CAS tokeniser
# has no quote handling); user name is "plug".
L.append('(45 (2 VC:porous plug)())')
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
# ---- y symmetry ----
f0, f1 = NF_INT + NF_ZSYM + 1, NF_INT + NF_ZSYM + NF_YSYM
L.append('(0 "Zone 5 %d faces %d..%d, symmetry")' % (NF_YSYM, f0, f1))
L.append('(13 (5 %s %s 7 4)(' % (hx(f0), hx(f1)))
for i in range(NX):   # y=0
    face([bn(i, 0), tn(i, 0), tn(i + 1, 0), bn(i + 1, 0)], cell(i, 0), 0)
for i in range(NX):   # y=NY
    face([bn(i, NY), tn(i, NY), tn(i + 1, NY), bn(i + 1, NY)], cell(i, NY - 1), 0)
L.append('))')
L.append('')
# ---- inlet x=0 ----
f0, f1 = NF_INT + NF_ZSYM + NF_YSYM + 1, NF_INT + NF_ZSYM + NF_YSYM + NF_IN
L.append('(0 "Zone 6 %d faces %d..%d, velocity-inlet")' % (NF_IN, f0, f1))
L.append('(13 (6 %s %s 5 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(0, j), tn(0, j), tn(0, j + 1), bn(0, j + 1)], cell(0, j), 0)
L.append('))')
L.append('')
# ---- outlet x=L ----
f0, f1 = NF - NF_OUT + 1, NF
L.append('(0 "Zone 7 %d faces %d..%d, pressure-outlet")' % (NF_OUT, f0, f1))
L.append('(13 (7 %s %s 4 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(NX, j), tn(NX, j), tn(NX, j + 1), bn(NX, j + 1)],
         cell(NX - 1, j), 0)
L.append('))')
L.append('')

with open(out, 'w') as f:
    f.write('\n'.join(L))
print('wrote', out, ': nodes', NN, 'cells', NC, 'faces', NF)
