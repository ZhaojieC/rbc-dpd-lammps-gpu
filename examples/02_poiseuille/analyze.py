"""Fit the velocity profile of in.poiseuille and report the viscosity.

For a body force f per particle, density rho and walls at z1 and z2,
  u(z) = rho f (z - z1) (z2 - z) / (2 eta).
"""

import numpy as np

F, RHO, Z1, Z2 = 0.5, 4.0, 1.5, 11.5

rows = []
for line in open("profile.txt"):
    w = line.split()
    if line.startswith("#") or len(w) != 5:
        continue
    rows.append([float(v) for v in w[1:]])
z, n, vx, dens = np.array(rows).T
core = (z > Z1 + 0.5) & (z < Z2 - 0.5)
c = np.polyfit(z[core], vx[core], 2)                     # free parabola
eta = -RHO * F / (2.0 * c[0])
roots = np.sort(np.roots(c).real)
print("viscosity eta = %.1f (reduced units)" % eta)
print("zero-velocity planes of the fit: z = %.2f and %.2f (walls at %.1f, %.1f)"
      % (roots[0], roots[1], Z1, Z2))
print("u_max = %.4f, mean density in the core %.3f" % (vx[core].max(), dens[core].mean()))
np.savetxt("fit.txt", np.c_[z, vx, np.polyval(c, z), dens],
           header="z vx fit density", fmt="%.5f")
