"""Speed, position and shape of the cell, from the output of in.channel.

  python3 analyze.py
"""

import numpy as np

VEL = 7.14                   # mm/s per DPD velocity unit (see the main README)


def last_frame(path):
    lines = open(path).read().splitlines()
    i = [k for k, l in enumerate(lines) if l.startswith("ITEM: ATOMS")][-1]
    step = int(lines[[k for k, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")][-1] + 1])
    return step, np.array([[float(v) for v in l.split()[2:5]] for l in lines[i + 1:]])


step, xc, zc, uc, area, vol = np.loadtxt("cell.txt", unpack=True)
late = step >= step[-1] / 2
print("cell speed (second half of the run) %.4f = %.2f mm/s" % (uc[late].mean(), uc[late].mean() * VEL))
print("height of the cell center: start %.2f, end %.2f (channel center 6.5)" % (zc[0], zc[-1]))
print("area change %+.2f %%, volume change %+.2f %%" % (100 * (area[-1] / area[0] - 1),
                                                     100 * (vol[-1] / vol[0] - 1)))
s, x = last_frame("cell.lammpstrj")
ext = x.max(axis=0) - x.min(axis=0)
print("cell size at step %d: %.2f (x, flow) x %.2f (y) x %.2f (z) um" % (s, *ext))
