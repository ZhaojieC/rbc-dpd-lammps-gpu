# 05_slit: a cell squeezing through a slit

A red cell pushed through a slit 1.2 um high and 5 um wide, cut in a
plate 2.5 um thick (box 24 x 12 x 13 um, periodic in x and y, 15,483
particles). Walls at the bottom and top leave a channel 10 um high, and
the plate closes it except for the slit. The cell starts flat, its
center 5.5 um upstream of the plate. A force on every membrane vertex
drives it through the slit; the solvent is not driven. The run stops
when the whole cell has passed the plate.

<p align="center"><img src="../../docs/slit.gif" width="640" alt="A red cell passing through a slit, side view"></p>

```
python3 make_data.py
../../build/lmp -k on g 1 -sf kk -pk kokkos newton on neigh half -in in.slit
python3 analyze.py
```

The run takes about 50,000 steps (7 ms). Run times are in the main README.

## How the cell is driven

In experiments a pressure difference drives cells through slits, and
the flow drags them along. The flow converges toward the slit, is
fastest inside it and slows down again past it: the same flux passes
through 120 um^2 of channel and through the 6 um^2 of the slit. The DPD
solvent here is compressible, so a body force on it would mostly pile
solvent up in front of the plate. Instead, a force on every membrane
vertex mimics that drag:

- upstream of the plate it points to the middle of the slit. Its size
  is f far from the plate and grows to `boost` x f over the last `w` um
  before the plate;
- inside the slit it is `boost` x f along x;
- past the slit it stays along x and falls back to f over `w2` um.

Here f = 40 (1.6 pN), `boost` = 3, `w` = 1 um and `w2` = 2 um. The
force is set by the atom-style variables `fx`, `fy` and `fz` in
`in.slit`. It is not calibrated to a pressure: use it to compare cells
or slits with each other, not to predict transit times in the spleen.

## The membrane slides on the walls

The membrane presses on the plate and on the walls of the slit, and it
has to slide along them to get through. Two settings let it slide.
`slip cell` in `walls.inc` (written by `make_data.py`) reflects the
membrane vertices like a mirror, so only their velocity normal to the
wall is reversed, and `pair_coeff 1 3` gives the membrane no drag from
the wall particles. The solvent keeps its no-slip walls, so a membrane
that slides along a wall still feels the drag of the still solvent
there. With no-slip walls for the membrane (no `slip cell`, and the same
drag from the wall particles as from the solvent), the cell crawls once
it is pressed on the plate: in our runs it moved 0.45 um from step 20,000
to step 30,000, against 2.0 um with slip.

The membrane also sees the edges of the slit rounded, with a radius of
0.3 um (`round 0.3` in `walls.inc`, `R` in `make_data.py`); the solvent
sees them sharp. The springs of the membrane are 0.4 to 0.7 um long. At
a sharp edge, a vertex can hang on the face of the plate just outside
the slit while the springs that run from it into the slit pull it
straight against the face. Those springs then stretch until LAMMPS stops
with `Bond ... reached its maximum length`. With sharp edges for the
membrane as well, one of our three runs with a 1.5 um plate stopped this
way. With the 2.5 um plate used here all three passed, but the longest
springs reached 2.06 to 2.07 times their rest length (the limit,
`lmax/l0` in `bond_coeff`, is 2.2), against 1.92 to 1.95 with rounded
edges. The passage also took longer: 59 to 60 DPD units instead of 47.

## What to expect

`analyze.py` printed, for one of our runs:

```
area   130.39 to 136.39 (start 132.87)
volume 89.03 to 92.45 (start 92.45)
front reaches the plate at step 3000 (t = 3, 0.4 ms)
rear leaves the plate at step 50000 (t = 50, 7.0 ms)
passage time 47 DPD units = 6.6 ms
```

All our runs gave a passage time of 47 DPD units (five random seeds,
on 4 or 12 MPI processes or on one GPU). `cell.txt` has a line every
1000 steps, so passage times are known only to 1 unit.

The area and the volume are kept by penalties (`ka` and `kv` in
`angle_coeff`), not exactly, so they give a little while the cell is
squeezed: here the area changed by -2 to +3 % and the volume dropped by
up to 4 %. A cell at rest (`01_rest`) keeps both within 0.2 %.

## Output

| File | Content |
|---|---|
| `cell.txt` | step, center, front and rear of the cell along x, area and volume, every 1000 steps |
| `cell.lammpstrj` | vertex positions, every 1000 steps |
| `slit.restart` | the final state |

## Try

- Change the force with `-var f`, `-var boost`, `-var w` and `-var w2`,
  for example `-var f 30`.
- Make the slit narrower or wider (`SW` and `SH` in `make_data.py`).
  The force points to the middle of the slit (`x1`, `x2`, `yc` and `zc`
  in `in.slit`); change them too if you move the slit or the plate, and
  `X1`, `X2` in `analyze.py`.
- Make the cell stiffer, for example with a larger shear modulus in
  `bond_coeff` (the second number), and compare the passage times.

If LAMMPS stops with `Bond ... reached its maximum length`, a spring of
the membrane was pulled past its limit: the push is too strong for the
slit. Lower `f` or `boost`, make the slit wider, or round its edges
more (a larger `R`).
