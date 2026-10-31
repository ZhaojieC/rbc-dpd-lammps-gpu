"""Shear rate, cell orientation and membrane rotation, from in.shear.

  python3 analyze.py

The orientation is the angle between the long axis of the cell in the
x-z plane and the flow direction x.  The membrane rotation rate is the
mean angular velocity of the vertices about the cell center in that
plane: a cell that tank-treads keeps its orientation while its membrane
rotates.
"""

import numpy as np

MS = 0.140                   # ms per DPD time unit (see the main README)


def profile(path):
    """Last block of a fix ave/chunk file: bin centers and vx."""
    rows, block = [], []
    for line in open(path):
        if line.startswith("#"):
            continue
        w = line.split()
        if len(w) == 3:
            block = []
            rows.append(block)
        else:
            block.append([float(v) for v in w])
    a = np.array(rows[-1])
    return a[:, 1], a[:, 3]


def frames(path):
    with open(path) as f:
        lines = f.read().splitlines()
    starts = [k for k, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    for k in starts:
        step = int(lines[k + 1])
        n = int(lines[k + 3])
        yield step, np.array([[float(v) for v in l.split()[2:5]] for l in lines[k + 9:k + 9 + n]])


z, vx = profile("profile.txt")
mid = (z > 4.5) & (z < 10.5)             # the middle half of the gap (1.5 < z < 13.5)
rate = np.polyfit(z[mid], vx[mid], 1)[0]
print("shear rate in the middle of the gap %.4f (%.0f 1/s)" % (rate, rate / MS * 1000))

steps, angle, phi = [], [], []
for step, x in frames("cell.lammpstrj"):
    r = x - x.mean(axis=0)
    g = np.cov(r[:, [0, 2]].T)
    w, v = np.linalg.eigh(g)
    e = v[:, 1]
    angle.append(np.degrees(np.arctan2(e[1], e[0])) % 180.0)
    phi.append(np.arctan2(r[:, 2], r[:, 0]))
    steps.append(step)
steps = np.array(steps)
angle = np.array(angle)
angle[angle > 90] -= 180
phi = np.unwrap(np.array(phi), axis=0)
t = steps * 0.001
late = steps >= steps[-1] / 2
omega = np.polyfit(t[late], phi[late].mean(axis=1), 1)[0]
print("orientation (second half): %.1f +- %.1f degrees" % (angle[late].mean(), angle[late].std()))
print("membrane rotation rate %.4f rad per time unit = %.3f x shear rate" % (-omega, -omega / rate))
