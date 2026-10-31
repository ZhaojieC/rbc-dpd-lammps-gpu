"""A red cell in pressure-driven flow between two walls.

Box 24 x 16 x 13 um, periodic in x and y; walls 1.5 um thick at the
bottom and top leave a channel 10 um high.  The cell starts at
mid-height with its disk facing the flow (axis along x).
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import rbcdata  # noqa: E402
from rbcdata import INF  # noqa: E402

LX, LY, H, W = 24.0, 16.0, 10.0, 1.5
walls = [(-INF, INF, -INF, INF, -INF, W),
         (-INF, INF, -INF, INF, W + H, INF)]

rbcdata.check_mesh()
s = rbcdata.System((0, LX, 0, LY, 0, H + 2 * W), seed=3)
s.add_cell(center=(6.0, LY / 2, W + H / 2), rot=rbcdata.rotation("y", 90))
s.add_particles(atype=2, density=4.0, exclude=walls)
for b in walls:
    s.add_particles(atype=3, density=4.0, region=b)
s.write("channel.data", "RBC in channel flow")
rbcdata.write_walls("walls.inc", walls, group="mobile")
