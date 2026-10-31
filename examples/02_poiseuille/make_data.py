"""DPD solvent between two flat walls (no cell).

Box 10 x 10 x 13 um: periodic in x and y; walls 1.5 um thick at the
bottom and top leave a channel 10 um high (1.5 < z < 11.5).
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import rbcdata  # noqa: E402
from rbcdata import INF  # noqa: E402

LX, LY, H, W = 10.0, 10.0, 10.0, 1.5
walls = [(-INF, INF, -INF, INF, -INF, W),
         (-INF, INF, -INF, INF, W + H, INF)]

s = rbcdata.System((0, LX, 0, LY, 0, H + 2 * W), seed=2)
s.add_particles(atype=2, density=4.0, exclude=walls)
for b in walls:
    s.add_particles(atype=3, density=4.0, region=b)
s.write("poiseuille.data", "DPD solvent between walls")
rbcdata.write_walls("walls.inc", walls, group="solvent")
