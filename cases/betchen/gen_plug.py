#!/usr/bin/env python3
# Betchen 2006 validation case 2 (PLUG): fluid / porous / fluid three-segment
# channel, interfaces PERPENDICULAR to the flow.
#
#   y in [0, H]            : channel half? No -- full height, H channel height.
#   x segments (default)   : dx1=3H (fluid) | dx2=2H (porous) | dx3=3H (fluid)
#   x segments (high-Re)   : dx1=5H | dx2=5H | dx3=50H  (long outlet section)
#
# Cell zones:
#   zone 2 fluid   x in [0, dx1]                     (VC:fluid)
#   zone 3 porous  x in [dx1, dx1+dx2]               (VC:porous)
#   zone 4 fluid   x in [dx1+dx2, dx1+dx2+dx3]       (VC:fluid)
# (A fluent cell-zone is one contiguous id range, so cid() must vary x SLOWEST.)
#
# Face zones: 5 interior, 6 z-symmetry, 7 wall y=0, 8 wall y=H,
#             9 velocity-inlet (x=0), 10 pressure-outlet (x=L).
import sys

H = 0.01                 # channel height (m)
if len(sys.argv) > 2 and sys.argv[2] == 'hir':
    dx1, dx2, dx3 = 5 * H, 5 * H, 50 * H     # high-Re variant
    NSEG = 20
else:
    dx1, dx2, dx3 = 3 * H, 2 * H, 3 * H       # default
    NSEG = 20
# per-segment cell counts in x (optional 3rd argument "n1,n2,n3")
NS1 = NS2 = NS3 = NSEG
if len(sys.argv) > 3 and sys.argv[3]:
    NS1, NS2, NS3 = [int(v) for v in sys.argv[3].split(',')]
LX = dx1 + dx2 + dx3
NX, NY = NS1 + NS2 + NS3, 21  # cells per segment in x, 21 in y (paper grid)
LZ = 0.002
out = sys.argv[1]
MX, MY = NX + 1, NY + 1
NP = MX * MY
NN = 2 * NP
rho, mu = 1.177, 1.846e-5
U0 = mu / (rho * H)      # Re_H = 1


def bn(i, j):
    return 1 + i + j * MX


def tn(i, j):
    return bn(i, j) + NP


def cell(i, j):
    # A Fluent cell-zone is one contiguous id range.  The 3 zones split by x,
    # so x must be the SLOWEST-varying cell index: all j for one i are
    # consecutive, then i advances.  (node ordering bn() keeps i fastest, but
    # cell ids are our own bookkeeping -- the face records reference them.)
    return 1 + j + i * NY


def xnode(i):
    if i < 0 or i > NX:
        raise ValueError
    # piecewise-uniform x: dx1/dx2/dx3 split into NS1/NS2/NS3 cells
    def xi(i0, L, n):
        return L * i0 / n
    if i <= NS1:
        return xi(i, dx1, NS1)
    if i <= NS1 + NS2:
        return dx1 + xi(i - NS1, dx2, NS2)
    return dx1 + dx2 + xi(i - NS1 - NS2, dx3, NS3)


def hx(v):
    return format(v, 'x')


nseg1, nseg2 = NS1, NS2            # cells in fluid-1, porous, fluid-3
nc1 = NS1 * NY                     # fluid zone 2
nc2 = NS2 * NY                     # porous zone 3
nc3 = NS3 * NY                     # fluid zone 4
NF_YINT = (NY - 1) * NX
NF_XINT = (NX - 1) * NY
NF_INT = NF_YINT + NF_XINT
NF_ZSYM = 2 * NX * NY
NF_W0 = NX
NF_WH = NX
NF_IN = NY
NF_OUT = NY
NF = NF_INT + NF_ZSYM + NF_W0 + NF_WH + NF_IN + NF_OUT
NC = NX * NY

L = []
L.append('(0 "generated Betchen PLUG channel (fluid/porous/fluid)")')
L.append('(0 "NX=%d NY=%d, seg=%gH/%gH/%gH, H=%.4f m")'
         % (NX, NY, dx1 / H, dx2 / H, dx3 / H, H))
L.append('(0 "Re_H=1 -> U0=%.6e m/s")' % U0)
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
    yj = H * j / NY
    for i in range(MX):
        L.append('  %.15e   %.15e   0.000000000000000e+00'
                 % (xnode(i), yj))
for j in range(MY):
    yj = H * j / NY
    for i in range(MX):
        L.append('  %.15e   %.15e   %.15e' % (xnode(i), yj, LZ))
L.append('))')
L.append('')
# three cell zones, x varies slowest -> contiguous id ranges
c2 = cell(0, 0)
z22 = cell(nseg1 - 1, NY - 1)
z30 = z22 + 1
z32 = cell(nseg1 + nseg2 - 1, NY - 1)
z40 = z32 + 1
z42 = cell(NX - 1, NY - 1)
L.append('(0 "Zone 2 %d cells %d..%d, fluid seg1")'
         % (nc1, c2, z22))
L.append('(12 (2 %s %s 1 4))' % (hx(c2), hx(z22)))
L.append('(45 (2 fluid seg1)())')
L.append('')
L.append('(0 "Zone 3 %d cells %d..%d, VC:porous seg2")'
         % (nc2, z30, z32))
L.append('(12 (3 %s %s 1 4))' % (hx(z30), hx(z32)))
L.append('(45 (3 VC:porous seg2)())')
L.append('')
L.append('(0 "Zone 4 %d cells %d..%d, fluid seg3")'
         % (nc3, z40, z42))
L.append('(12 (4 %s %s 1 4))' % (hx(z40), hx(z42)))
L.append('(45 (4 fluid seg3)())')
L.append('')


def face(ns, c0, c1):
    L.append(' '.join(hx(v) for v in ns) + ' ' + hx(c0) + ' ' + hx(c1))


f0, f1 = 1, NF_INT
L.append('(0 "Zone 5 %d faces %d..%d, Interior")' % (NF_INT, f0, f1))
L.append('(13 (5 %s %s 2 4)(' % (hx(f0), hx(f1)))
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
L.append('(0 "Zone 6 %d faces %d..%d, symmetry")' % (NF_ZSYM, f0, f1))
L.append('(13 (6 %s %s 7 4)(' % (hx(f0), hx(f1)))
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
f0, f1 = NF_INT + NF_ZSYM + 1, NF_INT + NF_ZSYM + NF_W0
L.append('(0 "Zone 7 %d faces %d..%d, wall y=0")' % (NF_W0, f0, f1))
L.append('(13 (7 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(NX):
    face([bn(i, 0), tn(i, 0), tn(i + 1, 0), bn(i + 1, 0)], cell(i, 0), 0)
L.append('))')
L.append('')
f0, f1 = NF_INT + NF_ZSYM + NF_W0 + 1, NF_INT + NF_ZSYM + NF_W0 + NF_WH
L.append('(0 "Zone 8 %d faces %d..%d, wall y=H")' % (NF_WH, f0, f1))
L.append('(13 (8 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(NX):
    face([bn(i, NY), tn(i, NY), tn(i + 1, NY), bn(i + 1, NY)],
         cell(i, NY - 1), 0)
L.append('))')
L.append('')
f0, f1 = NF_INT + NF_ZSYM + NF_W0 + NF_WH + 1, \
         NF_INT + NF_ZSYM + NF_W0 + NF_WH + NF_IN
L.append('(0 "Zone 9 %d faces %d..%d, velocity-inlet x=0")'
         % (NF_IN, f0, f1))
L.append('(13 (9 %s %s 5 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(0, j), tn(0, j), tn(0, j + 1), bn(0, j + 1)], cell(0, j), 0)
L.append('))')
L.append('')
f0, f1 = NF - NF_OUT + 1, NF
L.append('(0 "Zone 10 %d faces %d..%d, pressure-outlet x=L")'
         % (NF_OUT, f0, f1))
# NOTE: CAS (13 zone-id tokens are parsed AS HEX by mod_uns_cas_reader
# (tok_int uses (Z...)); single digits 5..9 are fine, but id "10" must be
# written "a" so it parses back to integer 10 (control file bc=10 is decimal).
L.append('(13 (a %s %s 4 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    face([bn(NX, j), tn(NX, j), tn(NX, j + 1), bn(NX, j + 1)],
         cell(NX - 1, j), 0)
L.append('))')
L.append('')

with open(out, 'w') as f:
    f.write('\n'.join(L))
print('wrote', out, ': nodes', NN, 'cells', NC, 'faces', NF,
      ' seg=%g/%g/%g H' % (dx1 / H, dx2 / H, dx3 / H))