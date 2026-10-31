"""Data for the force tests: one cell in a 20 um periodic box.

Writes cell.data (the reference shape) and, for each shift along x,
perturb_<shift>.inc, which moves every vertex by a random amount of up
to 0.05 um.  A shift of 9 puts the cell across the periodic boundary.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import rbcdata  # noqa: E402

rbcdata.check_mesh()
s = rbcdata.System((0, 20, 0, 20, 0, 20), seed=7)
x = s.add_cell(center=(10, 10, 10))
s.write("cell.data", "RBC force test")

dx = np.random.default_rng(11).uniform(-0.05, 0.05, size=x.shape)
for shift in (0, 9):
    with open("perturb_%d.inc" % shift, "w") as f:
        for i, (px, py, pz) in enumerate(x + dx):
            # wrapped position plus image flag, so the cell stays in one piece
            ix, wx = divmod(px + shift, 20.0)
            f.write("set atom %d x %.10f y %.10f z %.10f image %d 0 0\n" % (i + 1, wx, py, pz, ix))
