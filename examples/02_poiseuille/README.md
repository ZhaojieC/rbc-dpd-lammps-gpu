# 02_poiseuille: solvent viscosity

DPD solvent between two flat walls 10 um apart (box 10 x 10 x 13 um,
5,196 particles), with no cell. A body force f = 0.5 on every solvent
particle drives the flow along x. For a fluid of density rho and
viscosity eta between walls at z1 and z2, the velocity profile is the
parabola

    u(z) = rho f (z - z1) (z2 - z) / (2 eta)

and `analyze.py` fits it to the measured profile to get eta.

```
python3 make_data.py
../../build/lmp -k on g 1 -sf kk -pk kokkos newton on neigh half -in in.poiseuille
python3 analyze.py
```

The run is 20,000 steps to reach steady flow, then 30,000 steps of
averaging. Run times are in the main README.

## What to expect

`analyze.py` printed, for one of our runs:

```
viscosity eta = 219.7 (reduced units)
zero-velocity planes of the fit: z = 1.17 and 11.79 (walls at 1.5, 11.5)
u_max = 0.1280, mean density in the core 3.996
```

The viscosity matches the value of 218 used to set the time unit (1.2
mPa s, see the main README). The fitted profile goes to zero about
0.3 um beyond each wall: the walls are not quite no-slip. Use the
fitted planes, not the wall positions, when you need the effective
channel height.

## Output

| File | Content |
|---|---|
| `profile.txt` | velocity and density in 0.25 um slabs across the channel |
| `fit.txt` | written by `analyze.py`: z, measured and fitted velocity, density |
