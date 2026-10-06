#!/usr/bin/env python3
"""Split the merged symmetry zone of Pointwise channel.cas.

The exported file assigns BOTH the 2D extrusion planes (z=0, z=10) and the
physical duct walls (y=0, y=10) to a single Fluent BC zone 4 ("symmetry",
60600 faces).  A .control file can bind only one BC type per zone id, so the
wall faces cannot be no-slip while the extrusion planes stay slip.

Geometry audit of the zone-4 payload (boundary quads, implicit global face
ids assigned by record order; 30000 hex cells, 1-cell z extrusion):
    records  1..30000  -> faces  59601..89600  z=0 plane   (symmetry)
    records 30001..30600 -> faces 89601..90200 y=0/y=10    (wall, 300 each)
    records 30601..60600 -> faces 90201..120200 z=10 plane (symmetry)
The three slices are contiguous and already ordered this way, so the file
can be split into three (13 face sections with explicit first/last ranges;
no face is renumbered and the cell section (12) is touched.  The solver's
CAS reader stores faces at absolute slots (mesh_append_faces), and two
sections sharing zone id 4 are merged through zone_find.

Usage: python3 split_channel_zones.py [src.cas] [dst.cas]
Defaults: channel.cas -> channel_wall.cas in this directory.
"""
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'channel.cas')
DST = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'channel_wall.cas')

raw = open(SRC, 'rb').read()
lines = raw.splitlines(keepends=True)   # CRLF preserved

def find(pred, start=0):
    for i in range(start, len(lines)):
        if pred(lines[i]):
            return i
    raise RuntimeError('marker not found')

# ---- zone 4 section -------------------------------------------------------
h = find(lambda s: s.startswith(b'(13 (4 '))
assert lines[h].startswith(b'(13 (4 e8d1 1d588 7 4)('), lines[h]
e = h + 1
payload = []
while lines[e].strip() != b'))':
    payload.append(lines[e])
    e += 1
assert len(payload) == 60600, len(payload)
z0, yw, z10 = payload[:30000], payload[30000:30600], payload[30600:]

# sanity: count tokens == 6 on every record (boundary quad: v0..v3 c0 0)
for r in payload:
    assert len(r.split()) == 6, r

def hx(n):
    return ('%x' % n).encode()

def section(zid, first, last, cond, records):
    head = b'(13 (%s %s %s %d 4)(\r\n' % (hx(zid), hx(first), hx(last), cond)
    return [head] + records + [b'))\r\n']

# global face ids: z0 59601..89600, wall 89601..90200, z10 90201..120200
new_secs = []
new_secs += section(4, 59601, 89600, 7, z0)    # symmetry (cond 7)
new_secs += section(7, 89601, 90200, 3, yw)    # wall     (cond 3)
new_secs += section(4, 90201, 120200, 7, z10)  # symmetry

out = lines[:h] + new_secs + lines[e+1:]

# ---- descriptive comment lines (skipped by the parser, kept for humans) --
txt = b''.join(out)
txt = txt.replace(
    b'(0 "Zone 4 60600 faces 59601..120200, BC: bc-2 symmetry = 7")',
    b'(0 "Zone 4 60000 faces 59601..89600 + 90201..120200, BC: bc-2 symmetry = 7 (z planes, split)")\r\n'
    b'(0 "Zone 7 600 faces 89601..90200, BC: bc-wall wall = 3 (y=0/y=10, split from zone 4)")')

# ---- zone name record (45) for the new wall zone --------------------------
anchor = b'(45 (4 symmetry bc-2)())'
i = txt.find(anchor)
assert i > 0
end = txt.find(b'\r\n', i) + 2
txt = txt[:end] + b'(45 (7 wall bc-wall)())\r\n' + txt[end:]

open(DST, 'wb').write(txt)
print('wrote', DST, '(%.1f KB)' % (len(txt) / 1024))
