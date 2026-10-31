"""Passage of the cell through the slit, from cell.txt of in.slit.

  python3 analyze.py

The cell enters when its front reaches the upstream face of the plate
(x = 12) and has passed when its rear leaves the downstream face
(x = 14.5).  One DPD time unit is 0.14 ms (see the main README).
"""

import numpy as np

X1, X2 = 12.0, 14.5          # plate faces, as in make_data.py
DT = 0.001                   # time step
MS = 0.140                   # ms per DPD time unit

step, xc, front, rear, area, vol = np.loadtxt("cell.txt", unpack=True)
entered = step[front > X1]
passed = step[rear > X2]
print("area   %.2f to %.2f (start %.2f)" % (area.min(), area.max(), area[0]))
print("volume %.2f to %.2f (start %.2f)" % (vol.min(), vol.max(), vol[0]))
if len(entered) == 0:
    raise SystemExit("the cell has not reached the plate yet (front at x = %.2f)" % front[-1])
t_in = entered[0] * DT
print("front reaches the plate at step %d (t = %.0f, %.1f ms)" % (entered[0], t_in, t_in * MS))
if len(passed) == 0:
    raise SystemExit("the cell has not passed the plate yet (rear at x = %.2f)" % rear[-1])
t_out = passed[0] * DT
print("rear leaves the plate at step %d (t = %.0f, %.1f ms)" % (passed[0], t_out, t_out * MS))
print("passage time %.0f DPD units = %.1f ms" % (t_out - t_in, (t_out - t_in) * MS))
