"""Compare analytic and finite-difference forces from in.forces dumps.

  python compare.py forces_bond_0.dump [forces_bond_9.dump]

With one file: relative RMS difference between fx..fz and f_fd.
With two files: also the difference between the analytic forces of the
two runs (same cell, shifted across the periodic boundary).
"""

import sys

import numpy as np


def load(path):
    lines = open(path).read().splitlines()
    i = [k for k, l in enumerate(lines) if l.startswith("ITEM: ATOMS")][-1]
    return np.array([[float(v) for v in l.split()] for l in lines[i + 1:]])


def rel(a, b):
    return np.sqrt(((a - b) ** 2).sum() / (b ** 2).sum())


a = load(sys.argv[1])
err = rel(a[:, 1:4], a[:, 4:7])
print("%-26s |f| rms %10.4g   analytic vs numdiff %.2e" % (
    sys.argv[1], np.sqrt((a[:, 1:4] ** 2).sum(axis=1).mean()), err))
bad = err > 1e-6
if len(sys.argv) > 2:
    b = load(sys.argv[2])
    e2 = rel(b[:, 1:4], a[:, 1:4])
    print("%-26s shifted vs unshifted %.2e" % (sys.argv[2], e2))
    bad |= e2 > 1e-9
sys.exit(1 if bad else 0)
