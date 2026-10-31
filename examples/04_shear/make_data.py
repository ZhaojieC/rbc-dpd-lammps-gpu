"""A red cell in shear flow between two moving walls.

Box 16 x 8 x 15 um, periodic in x and y; walls 1.5 um thick at the
bottom and top leave a gap of 12 um.  The cell starts at the center
with its disk in the shear (x-z) plane.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import rbcdata  # noqa: E402
from rbcdata import INF  # noqa: E402

LX, LY, H, W = 16.0, 8.0, 12.0, 1.5
U = 0.6                                   # wall speed; shear rate 2 U / H = 0.1
walls = [(-INF, INF, -INF, INF, -INF, W),
         (-INF, INF, -INF, INF, W + H, INF)]

rbcdata.check_mesh()
s = rbcdata.System((0, LX, 0, LY, 0, H + 2 * W), seed=4)
s.add_cell(center=(LX / 2, LY / 2, W + H / 2), rot=rbcdata.rotation("x", 90))
s.add_particles(atype=2, density=4.0, exclude=walls)
for b in walls:
    s.add_particles(atype=3, density=4.0, region=b)
s.write("shear.data", "RBC in shear flow")
rbcdata.write_walls("walls.inc", walls, velocities=[(-U, 0, 0), (U, 0, 0)], group="mobile")
