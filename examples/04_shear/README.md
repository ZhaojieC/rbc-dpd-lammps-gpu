# 04_shear: a cell in shear flow

A red cell between two walls that slide in opposite directions (box
16 x 8 x 15 um, periodic in x and y, 8,178 particles). The walls are
12 um apart and move along x at -U and +U, with U = 0.6 (4.3 mm/s), so
the shear rate is 2U/12 = 0.1 (714 1/s). The wall particles stay in
place but carry the wall velocity, which their DPD drag passes on to the
solvent. The cell starts at the center as a disk lying in the plane of
the flow (x-z). It stretches and tilts, then keeps its shape while the
membrane goes round it like a tank tread.

<p align="center"><img src="../../docs/shear.gif" width="400" alt="A cell tank-treading in shear flow, viewed along the vorticity"></p>

The dark patch is a group of membrane vertices, colored only to show
the membrane moving.

```
python3 make_data.py
../../build/lmp -k on g 1 -sf kk -pk kokkos newton on neigh half -in in.shear
python3 analyze.py
```

The run is 250,000 steps (35 ms). Run times are in the main README.

## What to expect

`analyze.py` printed, for one of our runs:

```
shear rate in the middle of the gap 0.0761 (544 1/s)
orientation (second half): 13.3 +- 0.4 degrees
membrane rotation rate 0.0389 rad per time unit = 0.511 x shear rate
```

The cell keeps a tilt of about 13 degrees to the flow while its membrane
turns once every 160 time units (23 ms): it tank-treads. The shear rate
in the middle of the gap is lower than the 0.1 set by the walls, and
near the walls it is higher (about 0.11): the cell fills much of the gap,
and the flow is sheared less where the cell is.

## Output

| File | Content |
|---|---|
| `cell.txt` | step, center of the cell (x, z), area and volume, every 1000 steps |
| `profile.txt` | solvent velocity vx in 0.5 um layers across the gap, averaged over each 50,000 steps |
| `cell.lammpstrj` | vertex positions, every 2000 steps |
| `shear.restart` | the final state |

## Try

- Halve the shear rate: set `U = 0.3` in `make_data.py` and run it, then
  run LAMMPS with `-var U 0.3`. The wall speed is in both places, so
  change both.
- Start the cell with its disk across the flow (`rot=` in
  `make_data.py`) and follow its orientation in `cell.lammpstrj`.
