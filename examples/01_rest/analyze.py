"""Temperature, area, volume and size of the cell, from the output of in.rest.

  python3 analyze.py
"""

import numpy as np


def frames(path):
    with open(path) as f:
        lines = f.read().splitlines()
    for k, l in enumerate(lines):
        if l.startswith("ITEM: TIMESTEP"):
            n = int(lines[k + 3])
            yield int(lines[k + 1]), np.array([[float(v) for v in w.split()[2:5]]
                                               for w in lines[k + 9:k + 9 + n]])


step, temp, area, vol = np.loadtxt("shape.txt", unpack=True)
late = step >= step[-1] / 2
print("solvent temperature (second half) %.4f +- %.4f (set 0.0707)" % (temp[late].mean(), temp[late].std()))
print("area   %.2f at the start, largest change %.2f %%" % (area[0], 100 * np.abs(area / area[0] - 1).max()))
print("volume %.2f at the start, largest change %.2f %%" % (vol[0], 100 * np.abs(vol / vol[0] - 1).max()))
f = list(frames("cell.lammpstrj"))
for s, x in (f[0], f[-1]):
    print("size at step %6d: %.2f x %.2f x %.2f um" % (s, *(x.max(axis=0) - x.min(axis=0))))
