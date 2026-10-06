#!/usr/bin/env python3
# Generate the single-solver reference mesh for C2:
#
#   uns full    : Fluent ASCII .cas hex mesh, x in [0,200] mm
#                 two cell zones:
#                   zone 2 fluid      cells at x=0..100 mm  (i=0..49)
#                   zone 8 VC:porous  cells at x=100..200 mm (i=50..99)
#                 face zones:
#                   3 Interior, 4 symmetry (z), 5 symmetry (y),
#                   6 pressure-outlet (x=200), 7 velocity-inlet (x=0)
#
# This is the reference against which the coupled struct+uns run is judged,
# same methodology as couple_channel/uns_full.  All lateral faces are
# symmetry so the flow is strictly 1D.
#
# Units: millimetres; the .control sets mesh_scale=1e-3.
import struct

# ---- grid resolution (cells) -------------------------------------------------
NX, NY, NZ = 200, 40, 4
NF_HALF = 50                             # cells of fluid / porous each
LX, LY, LZ = 200.0, 50.0, 10.0           # mm, full domain
assert NF_HALF < NX

MX, MY, MZ = NX + 1, NY + 1, NZ + 1
NPNODE = MX * MY * MZ
NCELL = NX * NY * NZ

NF_X = (NX - 1) * NY * NZ
NF_Y = NX * (NY - 1) * NZ
NF_Z = NX * NY * (NZ - 1)
NF_INT = NF_X + NF_Y + NF_Z
NF_SYM_Z = 2 * NX * NY
NF_SYM_Y = 2 * NX * NZ
NF_OUT = NY * NZ
NF_IN = NY * NZ
NF = NF_INT + NF_SYM_Z + NF_SYM_Y + NF_OUT + NF_IN


def nid(i, j, k):          # 1-based node id
    return 1 + i + j * MX + k * MX * MY


def cid(i, j, k):          # 1-based cell id
    return 1 + i + j * NX + k * NX * NY


def hx(v):
    return format(v, 'x')


L = []
L.append('(0 "generated C2 single-solver reference (fluid + porous bed)")')
L.append('(0 "NX=%d NY=%d NZ=%d, mm units")' % (NX, NY, NZ))
L.append('')
L.append('(0 "Dimension : 3")')
L.append('(2 3)')
L.append('')
L.append('(0 "Number of Nodes : %d")' % NPNODE)
L.append('(10 (0 1 %s 0 3))' % hx(NPNODE))
L.append('')
L.append('(0 "Total Number of Faces : %d")' % NF)
L.append('(0 "       Boundary Faces : %d")'
         % (NF_SYM_Z + NF_SYM_Y + NF_OUT + NF_IN))
L.append('(0 "       Interior Faces : %d")' % NF_INT)
L.append('(13 (0 1 %s 0))' % hx(NF))
L.append('')
L.append('(0 "Total Number of Cells : %d")' % NCELL)
L.append('(0 "            Hex cells : %d")' % NCELL)
L.append('(12 (0 1 %s 0))' % hx(NCELL))
L.append('')

# ---- nodes ----
L.append('(0 "Zone 1  Number of Nodes : %d")' % NPNODE)
L.append('(10 (1 1 %s 1 3)(' % hx(NPNODE))
for k in range(MZ):
    for j in range(MY):
        for i in range(MX):
            L.append('  %.15e   %.15e   %.15e'
                     % (i * LX / NX, j * LY / NY, k * LZ / NZ))
L.append('))')
L.append('')

# ---- cells: two zones (fluid x<100 mm, VC:porous x>=100 mm) ------------------
NCZ = NF_HALF * NY * NZ
cz2_0, cz2_1 = 1, NCZ
cz8_0, cz8_1 = NCZ + 1, NCELL
L.append('(0 "Zone 2 %d cells %d..%d, fluid free stream")'
         % (NCZ, cz2_0, cz2_1))
L.append('(12 (2 %s %s 1 4))' % (hx(cz2_0), hx(cz2_1)))
L.append('(45 (2 fluid unspecified)())')
L.append('')
L.append('(0 "Zone 8 %d cells %d..%d, VC:porous bed")'
         % (NCZ, cz8_0, cz8_1))
L.append('(12 (8 %s %s 1 4))' % (hx(cz8_0), hx(cz8_1)))
L.append('(45 (8 VC:porous bed)())')
L.append('')


def face(ns, c0, c1):
    L.append(' '.join(hx(v) for v in ns) + ' ' + hx(c0) + ' ' + hx(c1))


# ---- zone 3: interior faces ---------------------------------------------------
f0, f1 = 1, NF_INT
L.append('(0 "Zone 3 %d faces %d..%d, Interior")' % (NF_INT, f0, f1))
L.append('(13 (3 %s %s 2 4)(' % (hx(f0), hx(f1)))
# x-normal faces between i-1 and i (includes the fluid/porous junction i=50)
for i in range(1, NX):
    for k in range(NZ):
        for j in range(NY):
            face([nid(i, j, k), nid(i, j + 1, k),
                  nid(i, j + 1, k + 1), nid(i, j, k + 1)],
                 cid(i - 1, j, k), cid(i, j, k))
# y-normal faces between j-1 and j
for j in range(1, NY):
    for k in range(NZ):
        for i in range(NX):
            face([nid(i, j, k), nid(i + 1, j, k),
                  nid(i + 1, j, k + 1), nid(i, j, k + 1)],
                 cid(i, j - 1, k), cid(i, j, k))
# z-normal faces between k-1 and k
for k in range(1, NZ):
    for j in range(NY):
        for i in range(NX):
            face([nid(i, j, k), nid(i + 1, j, k),
                  nid(i + 1, j + 1, k), nid(i, j + 1, k)],
                 cid(i, j, k - 1), cid(i, j, k))
L.append('))')
L.append('')

# ---- zone 4: symmetry z planes ------------------------------------------------
f0 = NF_INT + 1
f1 = NF_INT + NF_SYM_Z
L.append('(0 "Zone 4 %d faces %d..%d, symmetry")' % (NF_SYM_Z, f0, f1))
L.append('(13 (4 %s %s 7 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    for i in range(NX):
        face([nid(i, j, 0), nid(i + 1, j, 0),
              nid(i + 1, j + 1, 0), nid(i, j + 1, 0)], cid(i, j, 0), 0)
for j in range(NY):
    for i in range(NX):
        face([nid(i, j, NZ), nid(i + 1, j, NZ),
              nid(i + 1, j + 1, NZ), nid(i, j + 1, NZ)], cid(i, j, NZ - 1), 0)
L.append('))')
L.append('(45 (4 symmetry bc-2)())')
L.append('')

# ---- zone 5: symmetry y=0 and y=H ---------------------------------------------
f0 = f1 + 1
f1 = f0 + NF_SYM_Y - 1
L.append('(0 "Zone 5 %d faces %d..%d, symmetry")' % (NF_SYM_Y, f0, f1))
L.append('(13 (5 %s %s 7 4)(' % (hx(f0), hx(f1)))
for k in range(NZ):
    for i in range(NX):
        face([nid(i, 0, k), nid(i + 1, 0, k),
              nid(i + 1, 0, k + 1), nid(i, 0, k + 1)], cid(i, 0, k), 0)
for k in range(NZ):
    for i in range(NX):
        face([nid(i, NY, k), nid(i + 1, NY, k),
              nid(i + 1, NY, k + 1), nid(i, NY, k + 1)], cid(i, NY - 1, k), 0)
L.append('))')
L.append('(45 (5 symmetry bc-3)())')
L.append('')

# ---- zone 6: pressure outlet at x = 200 mm ------------------------------------
f0 = f1 + 1
f1 = f0 + NF_OUT - 1
L.append('(0 "Zone 6 %d faces %d..%d, BC: bc-4 pressure-outlet = 5")'
         % (NF_OUT, f0, f1))
L.append('(13 (6 %s %s 5 4)(' % (hx(f0), hx(f1)))
for k in range(NZ):
    for j in range(NY):
        face([nid(NX, j, k), nid(NX, j + 1, k),
              nid(NX, j + 1, k + 1), nid(NX, j, k + 1)],
             cid(NX - 1, j, k), 0)
L.append('))')
L.append('(45 (6 pressure-outlet bc-4)())')
L.append('')

# ---- zone 7: velocity inlet at x = 0 ------------------------------------------
f0 = f1 + 1
f1 = NF
L.append('(0 "Zone 7 %d faces %d..%d, BC: bc-5 velocity-inlet = 4")'
         % (NF_IN, f0, f1))
L.append('(13 (7 %s %s 4 4)(' % (hx(f0), hx(f1)))
for k in range(NZ):
    for j in range(NY):
        face([nid(0, j, k), nid(0, j + 1, k),
              nid(0, j + 1, k + 1), nid(0, j, k + 1)],
             cid(0, j, k), 0)
L.append('))')
L.append('(45 (7 velocity-inlet bc-5)())')
L.append('')

with open('unMesh.cas', 'w') as f:
    f.write('\n'.join(L))

print('uns full    : %d nodes, %d hex cells (zone2 fluid %d, zone8 porous %d)'
      % (NPNODE, NCELL, NCZ, NCZ))
