# 01_rest: a cell at rest

One red cell in a periodic box of DPD solvent, 16 x 16 x 10 um, with
10,740 particles and no flow. Start here: if the build and the
parameters are right, the cell keeps its shape, area and volume, and
the solvent keeps its temperature.

```
python3 make_data.py
../../build/lmp -k on g 1 -sf kk -pk kokkos newton on neigh half -in in.rest
python3 analyze.py
```

The run is 50,000 steps (7 ms). Run times are in the main README.

## What to expect

`analyze.py` printed, for one of our runs:

```
solvent temperature (second half) 0.0709 +- 0.0007 (set 0.0707)
area   132.87 at the start, largest change 0.13 %
volume 92.45 at the start, largest change 0.14 %
size at step      0: 7.79 x 7.79 x 2.58 um
size at step  50000: 7.98 x 7.85 x 2.52 um
```

Your numbers will differ a little (thermal noise), but the area and
volume should stay within a few tenths of a percent.

## Output

| File | Content |
|---|---|
| `shape.txt` | step, solvent temperature, area and volume of the cell, every 1000 steps |
| `cell.lammpstrj` | vertex positions, every 2000 steps |
| `rest.restart` | the final state; continue from it with `read_restart` |
