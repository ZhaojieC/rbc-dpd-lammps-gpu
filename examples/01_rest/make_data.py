"""One red cell in a periodic box of DPD solvent (16 x 16 x 10 um)."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import rbcdata  # noqa: E402

rbcdata.check_mesh()
s = rbcdata.System((0, 16, 0, 16, 0, 10), seed=1)
s.add_cell(center=(8, 8, 5))
s.add_particles(atype=2, density=4.0)
s.write("rest.data", "RBC at rest in DPD solvent")
