"""A red cell in front of a slit.

Box 24 x 12 x 13 um, periodic in x and y.  Walls 1.5 um thick at the
bottom and top leave a channel 10 um high.  A plate 2.5 um thick
(12 < x < 14.5) closes the channel except for an aperture 5 um wide
(y) and 1.2 um high (z) at mid-height.  The cell lies flat (axis along
z) with its center 5.5 um upstream of the plate.  The membrane slides
on the walls and sees their edges rounded (slip cell and round in
walls.inc; see README.md).
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import rbcdata  # noqa: E402
from rbcdata import INF  # noqa: E402

LX, LY, H, W = 24.0, 12.0, 10.0, 1.5      # box, channel height, wall thickness
X1, X2 = 12.0, 14.5                       # plate faces
SW, SH = 5.0, 1.2                         # aperture width (y) and height (z)
R = 0.3                                   # radius of the wall edges seen by the membrane
yc, zc = LY / 2, W + H / 2

walls = [
    (-INF, INF, -INF, INF, -INF, W),                         # bottom wall
    (-INF, INF, -INF, INF, W + H, INF),                      # top wall
    (X1, X2, -INF, INF, W, zc - SH / 2),                     # plate below the aperture
    (X1, X2, -INF, INF, zc + SH / 2, W + H),                 # plate above it
    (X1, X2, -INF, yc - SW / 2, zc - SH / 2, zc + SH / 2),   # plate beside it
    (X1, X2, yc + SW / 2, INF, zc - SH / 2, zc + SH / 2),
]

rbcdata.check_mesh()
s = rbcdata.System((0, LX, 0, LY, 0, H + 2 * W), seed=5)
s.add_cell(center=(X1 - 5.5, yc, zc))
s.add_particles(atype=2, density=4.0, exclude=walls)
for b in walls:
    s.add_particles(atype=3, density=4.0, region=b)
s.write("slit.data", "RBC in front of a slit")

# The blocks of walls.inc overlap where the plate meets the channel walls
# and beside the aperture, so that the membrane finds rounded edges only
# at the slit.
mirror = [
    (-INF, INF, -INF, INF, -INF, W),
    (-INF, INF, -INF, INF, W + H, INF),
    (X1, X2, -INF, INF, -INF, zc - SH / 2),
    (X1, X2, -INF, INF, zc + SH / 2, INF),
    (X1, X2, -INF, yc - SW / 2, -INF, INF),
    (X1, X2, yc + SW / 2, INF, -INF, INF),
]
rbcdata.write_walls("walls.inc", mirror, group="mobile", slip="cell", round=R)
