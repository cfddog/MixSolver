#!/usr/bin/env python3
"""
check_iface_alignment.py -- coupling_C_S interface pairing guard.

WHAT IS CHECKED
---------------
The weak-coupling driver pairs the structured and unstructured interface faces
by INDEX, not geometrically: mod_mix_driver's production path maps struct
Interface_List entry #c <-> uns exchange row #c (the geometric matcher
match_interfaces() is a test-only path).  The two orders are built
independently:

  * struct  (src/structured/mod_struct_grid.f90: register_bc_interfaces)
      loops gridgen generic:8 subfaces in bc3d.inp order, and inside one
      subface iterates `a` (first in-plane index) OUTER and `b` (second)
      INNER, i.e. for a j-face the order is i increasing (k fastest, here k
      has a single cell).
  * uns     (src/unstructured/mod_uns_bc.f90: build_bc interface group)
      walks m%f(:) in increasing face index, i.e. the unMesh.cas face order
      inside the interface zone.

If those two orders disagree, the coupling still runs and still looks
plausible, but every face exchanges with a permuted partner -- a silent,
hard-to-spot physics error.  This script therefore re-derives BOTH orders
from the two mesh files and fails loudly if face #c does not pair with
face #c.

Case-specific regression guards:
  * zone 4 (slab base) must be named "wall" -- the delivered cas tagged it
    "mass-flow-inlet" and our copy patches it (see README);
  * zone 7 must be the interface zone;
  * exactly one code-8 block in bc3d.inp and it must be a j-face (the struct
    registration order replicated below assumes face 2/5);
  * both interface surfaces must lie in the y=0 plane, carry 80 faces each
    and 0.1 m^2 total area.

USAGE
-----
    python3 check_iface_alignment.py [--case-dir DIR] [--tol-m 1e-6]
Defaults to the directory containing this script.  Exit code 0 = PASS.
"""

from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path

FAILS: list[str] = []
WARNS: list[str] = []


def fail(msg: str) -> None:
    FAILS.append(msg)
    print("FAIL: " + msg)


def warn(msg: str) -> None:
    WARNS.append(msg)
    print("WARN: " + msg)


def note(msg: str) -> None:
    print("  " + msg)


# ---------------------------------------------------------------------------
# Fortran unformatted (sequential) reader: 4-byte record length before and
# after every record, little endian (Mesh3d.x, Plot3D-style, Mesh_File_Format=0)
# ---------------------------------------------------------------------------
def fortran_records(path: Path) -> list:
    data = path.read_bytes()
    recs = []
    off = 0
    while off + 4 <= len(data):
        n = struct.unpack_from("<i", data, off)[0]
        if n < 0 or off + 8 + n > len(data):
            raise ValueError(f"{path}: bad record length {n} at offset {off}")
        trailer = struct.unpack_from("<i", data, off + 4 + n)[0]
        if trailer != n:
            raise ValueError(f"{path}: record markers disagree at offset {off}")
        recs.append(data[off + 4:off + 4 + n])
        off += 8 + n
    return recs


def parse_struct_mesh(path: Path):
    """Return (ni, nj, nk, node(i,j,k) -> (x, y, z)) in mesh units (mm)."""
    recs = fortran_records(path)
    if len(recs) != 3:
        raise ValueError(f"{path}: expected 3 records, found {len(recs)}")
    hdr = struct.unpack("<i", recs[0])[0]
    ni, nj, nk = struct.unpack("<3i", recs[1])
    n = ni * nj * nk
    if len(recs[2]) != 3 * n * 8:
        raise ValueError(f"{path}: coordinate record size mismatch")
    v = struct.unpack("<%dd" % (3 * n), recs[2])
    bx, by, bz = v[:n], v[n:2 * n], v[2 * n:]

    # identify the three blocks by their coordinate ranges (distinct for this
    # case: x 0..1000, y 0..500, z 0..100 mm).  This detects an X/Y/Z layout
    # change in the writer and catches a swapped mesh file.
    spans = [(max(b) - min(b), min(b), max(b)) for b in (bx, by, bz)]
    note(f"Mesh3d.x header word = {hdr}, nodes = {ni} x {nj} x {nk} = {n}")
    for nm, (s, lo, hi) in zip("XYZ", spans):
        note(f"  stored block {nm}: range [{lo:g}, {hi:g}] mm")
    got = tuple(round(spans[i][0], 6) for i in range(3))
    if got != (1000.0, 500.0, 100.0):
        fail(f"Mesh3d.x block ranges (x,y,z) = {got} mm, expected "
             f"(1000.0, 500.0, 100.0) -- wrong mesh file or layout change")

    def node(i: int, j: int, k: int):
        p = (i - 1) + (j - 1) * ni + (k - 1) * ni * nj
        return (bx[p], by[p], bz[p])

    # the struct interface must be the j=1 plane (y = 0)
    ys = {round(node(i, 1, 1)[1], 12) for i in range(1, ni + 1)}
    if len(ys) != 1 or abs(next(iter(ys))) > 1e-9:
        fail(f"struct j=1 plane is not the y=0 plane (y values {ys})")
    note(f"  struct j=1 plane: y = {next(iter(ys)):g} mm")
    return ni, nj, nk, node


# ---------------------------------------------------------------------------
# bc3d.inp (Gridgen generic format): every data line is
#     ib ie jb je kb ke bc
# ---------------------------------------------------------------------------
def parse_bc3d(path: Path):
    blocks = []
    for line in path.read_text().splitlines():
        t = line.split()
        if len(t) == 7 and all(x.lstrip("-").isdigit() for x in t):
            blocks.append(tuple(int(x) for x in t))
    return blocks


# ---------------------------------------------------------------------------
# unMesh.cas (Fluent / Pointwise ascii export)
# ---------------------------------------------------------------------------
def cas_ints_are_hex(text: str) -> bool:
    """Pointwise writes cas integers in hexadecimal.  Detect it from the
    (13 ...) zone records: if any integer token there carries a hex digit
    a-f, the whole file uses base 16 (cross-checked against the decimal
    (0 \"Zone ...\") comments in parse_cas)."""
    for line in text.splitlines():
        if line.startswith("(13 ("):
            toks = line.replace("(", " ").replace(")", " ").split()
            if any(any(c in "abcdefABCDEF" for c in t) for t in toks[2:]):
                return True
    return False


def parse_cas(path: Path):
    text = path.read_text()
    lines = text.splitlines()
    hex_ints = cas_ints_are_hex(text)
    note(f"unMesh.cas integer tokens: "
         f"{'hexadecimal' if hex_ints else 'decimal'}")

    def to_int(tok: str) -> int:
        return int(tok, 16) if hex_ints else int(tok)

    # ---- (0 "Zone ...") comments: decimal summary of every zone ----------
    zones = {}
    for raw in lines:
        line = raw.strip()
        if not line.startswith("(0 "):
            continue
        body = line[line.find('"') + 1:line.rfind('"')]
        if not body.startswith("Zone ") or ".." not in body:
            continue
        t = body.split()
        if len(t) < 5:
            continue
        # zone comments come in two shapes:
        #   "Zone 2 3200 cells 1..3200, VC: ..."      (volume zone)
        #   "Zone 1  Number of Nodes : 6642"          (node zone, no range)
        try:
            zid, nent = int(t[1]), int(t[2])
            first, last = (int(v.strip(" .,")) for v in t[4].split(".."))
        except ValueError:
            continue
        zones[zid] = dict(nent=nent, first=first, last=last, comment=body)

    # ---- (45 (zone cond user)()) records: zone names ---------------------
    names = {}
    for raw in lines:
        line = raw.strip()
        if line.startswith("(45 ("):
            t = line.replace("(", " ").replace(")", " ").split()
            names[int(t[1])] = t[2]

    # ---- (13 (zone first last type elemtype)( + face lines --------------
    faces = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("(13 (") and line.endswith("("):
            t = line.replace("(", " ").split()
            zid, first, last = to_int(t[1]), to_int(t[2]), to_int(t[3])
            nf = last - first + 1
            quads = []
            for k in range(i + 1, i + 1 + nf):
                quads.append([to_int(x) for x in lines[k].split()[:4]])
            faces[zid] = dict(first=first, last=last, quads=quads)
            z = zones.get(zid)
            if z and (z["nent"] != nf or z["first"] != first or
                      z["last"] != last):
                fail(f"zone {zid}: (13 ...) says {nf} faces {first}..{last} "
                     f"but the (0 \"Zone ...\") comment says {z['nent']} "
                     f"{z['first']}..{z['last']} -- hex/decimal detection is "
                     f"wrong")
            i += nf
        i += 1

    # ---- node coordinates: (10 (1 1 <n> 1 3)( then 3 floats per node -----
    nnode = None
    nodes = []
    for i, raw in enumerate(lines):
        line = raw.strip()
        if line.startswith("(10 (1 ") and line.endswith("("):
            # token layout: (10 (1 1 <nhex> 1 3)( -> ['10','1','1','<nhex>',...]
            nnode = to_int(line.replace("(", " ").split()[3])
            for k in range(i + 1, i + 1 + nnode):
                x, y, z = (float(v) for v in lines[k].split()[:3])
                nodes.append((x, y, z))
            break
    if nnode is None or len(nodes) != nnode:
        raise ValueError(f"{path}: node block not found / wrong size")
    note(f"  nodes = {nnode}; face zones = {sorted(faces)}")
    return dict(hex_ints=hex_ints, zones=zones, names=names,
                faces=faces, nodes=nodes)


# ---------------------------------------------------------------------------
def quad_centroid_area_normal(p0, p1, p2, p3):
    """Planar quad -> (centroid, area, unit normal).  Triangulated as
    (P0,P1,P2)+(P0,P2,P3) so a slightly warped quad still gets a sane area."""
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    def norm(a):
        return (a[0] ** 2 + a[1] ** 2 + a[2] ** 2) ** 0.5

    c = [0.0, 0.0, 0.0]
    area = 0.0
    for (a, b, d) in ((p0, p1, p2), (p0, p2, p3)):
        n = cross(sub(b, a), sub(d, a))
        tri = 0.5 * norm(n)
        area += tri
        for i in range(3):
            c[i] += tri * (a[i] + b[i] + d[i]) / 3.0
    if area > 0.0:
        c = [v / area for v in c]
    n = cross(sub(p1, p0), sub(p2, p0))
    ln = norm(n)
    n = [v / ln for v in n] if ln > 0 else [0.0, 0.0, 0.0]
    return tuple(c), area, tuple(n)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="coupling_C_S interface alignment / pairing guard")
    ap.add_argument("--case-dir", default=str(Path(__file__).resolve().parent))
    ap.add_argument("--tol-m", type=float, default=1e-6,
                    help="max allowed per-face centroid mismatch [m]")
    args = ap.parse_args()
    d = Path(args.case_dir)

    print("=" * 74)
    print("coupling_C_S interface alignment check")
    print(f"  case dir : {d}")
    print("=" * 74)

    # ---------------- structured side ----------------
    print("--- structured side (Mesh3d.x + bc3d.inp) ---")
    ni, nj, nk, node = parse_struct_mesh(d / "Mesh3d.x")
    blocks = parse_bc3d(d / "bc3d.inp")
    itf = [b for b in blocks if b[6] == 8]
    note(f"bc3d.inp bc blocks: {len(blocks)}, code-8 (generic) blocks: "
         f"{len(itf)}")
    if len(itf) != 1:
        fail(f"expected exactly one code-8 block in bc3d.inp, found {len(itf)}")
        return 1

    def face_of(ib, ie, jb, je, kb, ke):
        if sum(1 for a, b in ((ib, ie), (jb, je), (kb, ke)) if a == b) != 1:
            return None
        if ib == ie:
            return 1 if ib == 1 else 4
        if jb == je:
            return 2 if jb == 1 else 5
        return 3 if kb == 1 else 6

    ib, ie, jb, je, kb, ke, _ = itf[0]
    fa = face_of(ib, ie, jb, je, kb, ke)
    fname = {1: "i-", 4: "i+", 2: "j-", 5: "j+", 3: "k-", 6: "k+"}.get(fa, "?")
    note(f"code-8 block: i {ib}..{ie}, j {jb}..{je}, k {kb}..{ke} -> face {fa}"
         f" ({fname} plane)")
    if fa not in (2, 5):
        fail("this script replicates the struct registration order for j-faces "
             "only; re-derive it for the new face index before trusting it")
        return 1

    # replicate register_bc_interfaces (mod_struct_grid.f90:1046-1069):
    # a = first in-plane index OUTER, b = second INNER
    if fa in (1, 4):
        i0, a0, a1, b0, b1, vary = ib, jb, je, kb, ke, "jk"
    else:
        i0, a0, a1, b0, b1, vary = jb, ib, ie, kb, ke, "ik"
    struct_faces = []
    for a in range(a0, a1):
        for b in range(b0, b1):
            if vary == "jk":
                cc = [(i0, a, b), (i0, a + 1, b), (i0, a + 1, b + 1),
                      (i0, a, b + 1)]
            else:
                cc = [(a, i0, b), (a + 1, i0, b), (a + 1, i0, b + 1),
                      (a, i0, b + 1)]
            struct_faces.append([node(*v) for v in cc])
    note(f"struct interface quads (registration order) = {len(struct_faces)}")
    if len(struct_faces) != 80:
        fail(f"expected 80 struct interface faces, got {len(struct_faces)}")

    # ---------------- unstructured side ----------------
    print("--- unstructured side (unMesh.cas) ---")
    cas = parse_cas(d / "unMesh.cas")
    if cas["names"].get(7) != "interface":
        fail(f"cas zone 7 is named '{cas['names'].get(7)}', expected "
             f"'interface' (the coupling exchange zone)")
    if cas["names"].get(4) != "wall":
        fail(f"cas zone 4 (slab base) is named '{cas['names'].get(4)}' -- the "
             f"delivered case tagged it 'mass-flow-inlet'; coupling_C_S patches "
             f"it to 'wall' (see README)")
    if 7 not in cas["faces"]:
        fail("cas has no zone-7 face list")
        return 1
    uns_faces = [[cas["nodes"][nid - 1] for nid in q]
                 for q in cas["faces"][7]["quads"]]
    note(f"uns interface quads (cas face order) = {len(uns_faces)} "
         f"(faces {cas['faces'][7]['first']}..{cas['faces'][7]['last']})")
    if len(uns_faces) != 80:
        fail(f"expected 80 uns interface faces, got {len(uns_faces)}")

    # ---------------- face-by-face pairing ----------------
    print("--- face-by-face pairing (struct registration order vs cas order) ---")
    S = [quad_centroid_area_normal(*f) for f in struct_faces]
    U = [quad_centroid_area_normal(*f) for f in uns_faces]
    n = min(len(S), len(U))
    worst, worst_idx = 0.0, -1
    area_s = area_u = 0.0
    dot_min = 1.0
    for k in range(n):
        dm = max(abs(S[k][0][i] - U[k][0][i]) for i in range(3)) / 1000.0
        if dm > worst:
            worst, worst_idx = dm, k + 1
        area_s += S[k][1] / 1e6
        area_u += U[k][1] / 1e6
        dot_min = min(dot_min, sum(S[k][2][i] * U[k][2][i] for i in range(3)))
    if worst > args.tol_m:
        fail(f"worst per-face centroid mismatch = {worst:.3e} m at face "
             f"{worst_idx} (> tol {args.tol_m:g}): the struct and uns face "
             f"orders do NOT coincide -- the coupled run would exchange "
             f"permuted partners, do not trust it")
    elif worst > 1e-9:
        warn(f"worst per-face centroid mismatch = {worst:.3e} m: index pairing "
             f"holds, but the meshes are not node-identical")
    note(f"worst per-face centroid mismatch : {worst:.3e} m (face {worst_idx})")
    note(f"sum area struct = {area_s:.10f} m^2, uns = {area_u:.10f} m^2")
    if abs(area_s - 0.1) > 1e-9:
        fail(f"struct interface area {area_s:.10f} m^2 != 0.1 m^2")
    if abs(area_u - area_s) > 1e-9:
        fail(f"interface areas differ: struct {area_s:.10f}, uns {area_u:.10f}")
    note(f"min normal dot product (struct vs uns) = {dot_min:+.6f}")
    if dot_min < 0.0:
        warn("some normals point the opposite way between the two meshes "
             "(harmless for the pairing, but watch the flux sign)")
    for idx in (0, 1, n // 2, n - 1):
        note(f"face {idx + 1:3d}: struct "
             f"{tuple(round(v, 9) for v in S[idx][0])} / uns "
             f"{tuple(round(v, 9) for v in U[idx][0])} (mm)")

    # ---------------- near-interface resolution ----------------
    print("--- near-interface resolution (CHT gain inputs) ---")
    dy_f = [node(1, j + 1, 1)[1] - node(1, j, 1)[1] for j in range(1, 6)]
    note("struct dy from the interface inward [mm]: "
         + ", ".join(f"{v:.6g}" for v in dy_f))
    x0, y0, z0 = uns_faces[0][0]
    col = sorted((p[1] for p in cas["nodes"]
                  if abs(p[0] - x0) < 1e-9 and abs(p[2] - z0) < 1e-9),
                 reverse=True)
    if len(col) >= 3:
        dy_s = [round(col[i] - col[i + 1], 9) for i in range(3)]
        note("uns slab dy from the interface inward [mm]: "
             + ", ".join(f"{v:.6g}" for v in dy_s))
        note(f"uns slab thickness = {col[0] - col[-1]:g} mm, "
             f"first layer = {dy_s[0] * 1000:.3f} um")
    note(f"struct first layer = {dy_f[0] * 1000:.3f} um")

    # ---------------- verdict ----------------
    print("=" * 74)
    if FAILS:
        print(f"RESULT: FAIL ({len(FAILS)} failure(s), {len(WARNS)} warning(s))")
        for f in FAILS:
            print("  - " + f)
        return 1
    print(f"RESULT: PASS (0 failures, {len(WARNS)} warning(s))")
    print("  the two interface face orders coincide 1:1 and the assets are the "
          "expected ones")
    return 0


if __name__ == "__main__":
    sys.exit(main())
