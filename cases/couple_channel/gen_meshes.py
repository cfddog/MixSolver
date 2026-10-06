#!/usr/bin/env python3
# Generate the coupled low-speed channel case meshes:
#
#   struct side : OpenCFD Plot3D unformatted Mesh3d.x  (x in [0,100] mm)
#                 + Gridgen bc3d.inp (i- forced inlet, i+ generic:8 interface,
#                 y walls, z symmetry)
#   uns side    : Fluent ASCII .cas hex mesh          (x in [100,200] mm)
#                 zone 4 symmetry (z), 5 wall (y),
#                 6 pressure-outlet (x=200), 7 interface (x=100)
#
# The two domains butt at x=100 mm.  The interface patch is identical on
# both sides (same 20x2 node grid, same face ORDER: j outer, k inner) so the
# phase-6 exchange, which maps faces by natural order, pairs them 1:1.
#
# Units: millimetres on both sides.  The unstructured .control sets
# mesh_scale=1e-3 so the uns solver internally works in metres.
import struct
import sys

# ---- grid resolution (cells) -------------------------------------------------
NX, NY, NZ = 50, 20, 2
LX, LY, LZ = 100.0, 50.0, 10.0          # mm, per domain
X0_UNS = 100.0                           # uns domain starts at x=100 mm

# ---- 1. structured side: Plot3D unformatted multi-block ----------------------
# Fortran sequential unformatted: each record is  [i32 len][payload][i32 len]
# with gfortran default 4-byte markers, little-endian.
NI, NJ, NK = NX + 1, NY + 1, NZ + 1

xs = [(i * LX / NX) for i in range(NI)]
ys = [(j * LY / NY) for j in range(NJ)]
zs = [(k * LZ / NZ) for k in range(NK)]


def frec(payload):
    return struct.pack('<i', len(payload)) + payload + struct.pack('<i', len(payload))


with open('Mesh3d.x', 'wb') as f:
    f.write(frec(struct.pack('<i', 1)))
    f.write(frec(struct.pack('<3i', NI, NJ, NK)))
    payload = bytearray()
    # x component, then y, then z; Fortran order (i fastest, j, k)
    for k in range(NK):
        for j in range(NJ):
            for i in range(NI):
                payload += struct.pack('<d', xs[i])
    for k in range(NK):
        for j in range(NJ):
            for i in range(NI):
                payload += struct.pack('<d', ys[j])
    for k in range(NK):
        for j in range(NJ):
            for i in range(NI):
                payload += struct.pack('<d', zs[k])
    f.write(frec(bytes(payload)))

# ---- bc3d.inp (single block, 6 physical faces, no internal connections) ------
# face codes: 2 wall, 3 symmetry, 5 inflow (forced), 8 coupling interface
bc = []
bc.append('     1')                              # ignored line
bc.append('     1')                              # number of blocks
bc.append('  %5d%5d%5d' % (NI, NJ, NK))
bc.append('channel-struct')
bc.append('     6')                              # 6 sub-faces
# i-  forced inflow
bc.append('  %6d%6d%6d%6d%6d%6d%6d' % (1, 1, 1, NJ, 1, NK, 5))
# i+  coupling interface
bc.append('  %6d%6d%6d%6d%6d%6d%6d' % (NI, NI, 1, NJ, 1, NK, 8))
# j-  wall
bc.append('  %6d%6d%6d%6d%6d%6d%6d' % (1, NI, 1, 1, 1, NK, 2))
# j+  wall
bc.append('  %6d%6d%6d%6d%6d%6d%6d' % (1, NI, NJ, NJ, 1, NK, 2))
# k-  symmetry
bc.append('  %6d%6d%6d%6d%6d%6d%6d' % (1, NI, 1, NJ, 1, 1, 3))
# k+  symmetry
bc.append('  %6d%6d%6d%6d%6d%6d%6d' % (1, NI, 1, NJ, NK, NK, 3))
with open('bc3d.inp', 'w') as f:
    f.write('\n'.join(bc) + '\n')

# ---- 2. unstructured side: Fluent ASCII CAS (structured hex grid) ------------
MX, MY, MZ = NX + 1, NY + 1, NZ + 1
NPNODE = MX * MY * MZ
NCELL = NX * NY * NZ

NF_X = (NX - 1) * NY * NZ
NF_Y = NX * (NY - 1) * NZ
NF_Z = NX * NY * (NZ - 1)
NF_INT = NF_X + NF_Y + NF_Z
NF_SYM = 2 * NX * NY
NF_WAL = 2 * NX * NZ
NF_OUT = NY * NZ
NF_IFC = NY * NZ
NF = NF_INT + NF_SYM + NF_WAL + NF_OUT + NF_IFC


def nid(i, j, k):          # 1-based node id
    return 1 + i + j * MX + k * MX * MY


def cid(i, j, k):          # 1-based cell id
    # NB: x is the FASTEST index here.  Harmless for THIS mesh because the uns
    # side is a SINGLE cell zone ("Zone 2", see unMesh.control: cell_zone = 2
    # fluid): the whole domain shares one zone, so no physics depends on the
    # ordering.  But a Fluent cell-zone record is one CONTIGUOUS id range: if
    # this mesh ever grows a second zone (an x-split fluid/porous pair), x must
    # become the SLOWEST index -- see the bugfix note in
    # couple_porous/uns_full/gen_full.py:cid(), where the k-slowest form
    # silently turned an intended x-SERIES zone split into two parallel z-slabs
    # (~1.8 Pa/mm instead of the analytic 3.02 Pa/mm).
    return 1 + i + j * NX + k * NX * NY


def hx(v):
    return format(v, 'x')


L = []
L.append('(0 "generated coupled low-speed channel (uns domain)")')
L.append('(0 "NX=%d NY=%d NZ=%d, mm units")' % (NX, NY, NZ))
L.append('')
L.append('(0 "Dimension : 3")')
L.append('(2 3)')
L.append('')
L.append('(0 "Number of Nodes : %d")' % NPNODE)
L.append('(10 (0 1 %s 0 3))' % hx(NPNODE))
L.append('')
L.append('(0 "Total Number of Faces : %d")' % NF)
L.append('(0 "       Boundary Faces : %d")' % (NF_SYM + NF_WAL + NF_OUT + NF_IFC))
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
            x = X0_UNS + i * LX / NX
            L.append('  %.15e   %.15e   %.15e' % (x, j * LY / NY, k * LZ / NZ))
L.append('))')
L.append('')

# ---- cells (no payload: connectivity reconstructed from faces) ----
L.append('(0 "Zone 2 %d cells 1..%d, fluid")' % (NCELL, NCELL))
L.append('(12 (2 1 %s 1 4))' % hx(NCELL))
L.append('(45 (2 fluid unspecified)())')
L.append('')


def face(ns, c0, c1):
    L.append(' '.join(hx(v) for v in ns) + ' ' + hx(c0) + ' ' + hx(c1))


# ---- zone 3: interior faces ---------------------------------------------------
f0, f1 = 1, NF_INT
L.append('(0 "Zone 3 %d faces %d..%d, Interior")' % (NF_INT, f0, f1))
L.append('(13 (3 %s %s 2 4)(' % (hx(f0), hx(f1)))
# x-normal faces between i-1 and i
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
f1 = NF_INT + NF_SYM
L.append('(0 "Zone 4 %d faces %d..%d, symmetry")' % (NF_SYM, f0, f1))
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

# ---- zone 5: walls y=0 and y=H ------------------------------------------------
f0 = f1 + 1
f1 = f0 + NF_WAL - 1
L.append('(0 "Zone 5 %d faces %d..%d, wall")' % (NF_WAL, f0, f1))
L.append('(13 (5 %s %s 3 4)(' % (hx(f0), hx(f1)))
for k in range(NZ):
    for i in range(NX):
        face([nid(i, 0, k), nid(i + 1, 0, k),
              nid(i + 1, 0, k + 1), nid(i, 0, k + 1)], cid(i, 0, k), 0)
for k in range(NZ):
    for i in range(NX):
        face([nid(i, NY, k), nid(i + 1, NY, k),
              nid(i + 1, NY, k + 1), nid(i, NY, k + 1)], cid(i, NY - 1, k), 0)
L.append('))')
L.append('(45 (5 wall bc-3)())')
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

# ---- zone 7: coupling interface at x = 100 mm ---------------------------------
# face order MUST match the struct registration: j outer, k inner.
f0 = f1 + 1
f1 = NF
L.append('(0 "Zone 7 %d faces %d..%d, BC: bc-5 interface = 24")'
         % (NF_IFC, f0, f1))
L.append('(13 (7 %s %s 24 4)(' % (hx(f0), hx(f1)))
for j in range(NY):
    for k in range(NZ):
        face([nid(0, j, k), nid(0, j + 1, k),
              nid(0, j + 1, k + 1), nid(0, j, k + 1)],
             cid(0, j, k), 0)
L.append('))')
L.append('(45 (7 interface bc-5)())')
L.append('')

with open('unMesh.cas', 'w') as f:
    f.write('\n'.join(L))

print('struct mesh : %dx%dx%d nodes, interface = %d quads' % (NI, NJ, NK, NY * NZ))
print('uns mesh    : %d nodes, %d hex cells, %d faces (interface %d)'
      % (NPNODE, NCELL, NF, NF_IFC))
