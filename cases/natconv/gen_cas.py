#!/usr/bin/env python3
# Generate a structured N x N x 1 hex .cas mesh identical in layout to
# cavity_natural_convection_porous.cas: 1m cube, 1 wedge-thick layer,
# zone2 fluid cells, zone3 interior faces, zone4 symmetry (z planes),
# zone5 wall (x=0,x=1,y=0,y=1). All indices hexadecimal.
import sys

N = int(sys.argv[1])
out = sys.argv[2]
M = N + 1            # nodes per side on one plane
NP = M * M           # nodes per plane
NN = 2 * NP          # total nodes

def bn(i, j):        # bottom-plane node id (z=0), 1-based
    return 1 + i + j * M

def tn(i, j):        # top-plane node id (z=1)
    return bn(i, j) + NP

def cell(i, j):
    return 1 + i + j * N

def hx(v):
    return format(v, 'x')

lines = []
lines.append('(0 "generated structured cavity")')
lines.append('(0 "N=%d x %d x 1")' % (N, N))
lines.append('')
lines.append('(0 "Dimension : 3")')
lines.append('(2 3)')
lines.append('')
lines.append('(0 "Number of Nodes : %d")' % NN)
lines.append('(10 (0 1 %s 0 3))' % hx(NN))
lines.append('')
NF_INT = 2 * N * (N - 1)
NF_SYM = 2 * N * N
NF_WAL = 4 * N
NF = NF_INT + NF_SYM + NF_WAL
lines.append('(0 "Total Number of Faces : %d")' % NF)
lines.append('(0 "       Boundary Faces : %d")' % (NF_SYM + NF_WAL))
lines.append('(0 "       Interior Faces : %d")' % NF_INT)
lines.append('(13 (0 1 %s 0))' % hx(NF))
lines.append('')
NC = N * N
lines.append('(0 "Total Number of Cells : %d")' % NC)
lines.append('(0 "            Hex cells : %d")' % NC)
lines.append('(12 (0 1 %s 0))' % hx(NC))
lines.append('')
# ---- nodes ----
lines.append('(0 "Zone 1  Number of Nodes : %d")' % NN)
lines.append('(10 (1 1 %s 1 3)(' % hx(NN))
for j in range(M):
    for i in range(M):
        lines.append('  %.15e   %.15e   0.000000000000000e+00' % (i / N, j / N))
for j in range(M):
    for i in range(M):
        lines.append('  %.15e   %.15e   1.000000000000000e+00' % (i / N, j / N))
lines.append('))')
lines.append('')
# ---- cells (uniform hex, no payload) ----
lines.append('(0 "Zone 2 %d cells 1..%d, fluid")' % (NC, NC))
lines.append('(12 (2 1 %s 1 4))' % hx(NC))
lines.append('(45 (2 fluid unspecified)())')
lines.append('')

def face(ns, c0, c1):
    lines.append(' '.join(hx(v) for v in ns) + ' ' + hx(c0) + ' ' + hx(c1))

# ---- interior faces ----
f0, f1 = 1, NF_INT
lines.append('(0 "Zone 3 %d faces %d..%d, Interior")' % (NF_INT, f0, f1))
lines.append('(13 (3 %s %s 2 4)(' % (hx(f0), hx(f1)))
# faces with normal in y: between rows j-1 and j
for j in range(1, N):
    for i in range(N):
        face([bn(i, j), tn(i, j), tn(i + 1, j), bn(i + 1, j)],
             cell(i, j - 1), cell(i, j))
# faces with normal in x: between columns i-1 and i
for i in range(1, N):
    for j in range(N):
        face([bn(i, j), tn(i, j), tn(i, j + 1), bn(i, j + 1)],
             cell(i - 1, j), cell(i, j))
lines.append('))')
lines.append('')
# ---- symmetry faces: z=0 and z=1 planes ----
f0 = NF_INT + 1
f1 = NF_INT + NF_SYM
lines.append('(0 "Zone 4 %d faces %d..%d, symmetry")' % (NF_SYM, f0, f1))
lines.append('(13 (4 %s %s 7 4)(' % (hx(f0), hx(f1)))
for j in range(N):
    for i in range(N):
        face([bn(i, j), bn(i + 1, j), bn(i + 1, j + 1), bn(i, j + 1)],
             cell(i, j), 0)
for j in range(N):
    for i in range(N):
        face([tn(i, j), tn(i + 1, j), tn(i + 1, j + 1), tn(i, j + 1)],
             cell(i, j), 0)
lines.append('))')
lines.append('')
# ---- wall faces: y=0,y=N,x=0,x=N ----
f0 = NF_INT + NF_SYM + 1
f1 = NF
lines.append('(0 "Zone 5 %d faces %d..%d, wall")' % (NF_WAL, f0, f1))
lines.append('(13 (5 %s %s 3 4)(' % (hx(f0), hx(f1)))
for i in range(N):   # y=0 hot
    face([bn(i, 0), tn(i, 0), tn(i + 1, 0), bn(i + 1, 0)], cell(i, 0), 0)
for i in range(N):   # y=N cold
    face([bn(i, N), tn(i, N), tn(i + 1, N), bn(i + 1, N)], cell(i, N - 1), 0)
for j in range(N):   # x=0
    face([bn(0, j), tn(0, j), tn(0, j + 1), bn(0, j + 1)], cell(0, j), 0)
for j in range(N):   # x=N
    face([bn(N, j), tn(N, j), tn(N, j + 1), bn(N, j + 1)], cell(N - 1, j), 0)
lines.append('))')
lines.append('')

with open(out, 'w') as f:
    f.write('\n'.join(lines))
print('wrote', out, ': nodes', NN, 'cells', NC, 'faces', NF)
