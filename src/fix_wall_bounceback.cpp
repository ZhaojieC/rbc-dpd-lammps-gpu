/* ----------------------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.

   Bounce-back (no-slip) reflection of particles at solid walls made of
   axis-aligned blocks.  Replaces fix bc/bbarb of the in-house code
   lmp_rbc (contributing author Igor Pivkin).

   After the position update of velocity Verlet, a particle found inside
   a block is moved back along its path: it travels to the wall, reverses
   its velocity relative to the wall (v -> 2 u_wall - v) and covers the
   rest of the step with the new velocity.  The old position x - dt v is
   exact for fix nve.  If the new position is still inside a block (near
   corners), the particle stays at its old position.

   Particles of the slip group are reflected like a mirror instead: only
   the velocity normal to the face they crossed is reversed, so they can
   slide along the wall.  Membrane vertices need this; with the full
   reversal, vertices pressed against a wall cannot move along it.

   With round R, the slip group sees the edges and corners of every block
   rounded with radius R (the block becomes the points closer than R to
   the block shrunk by R) and is mirrored at that surface.  A coarse
   membrane mesh then slides around the edge of a slit instead of
   hanging on it.

   fix ID group wall/bounceback [slip group-ID [round R]]
                                block xlo xhi ylo yhi zlo zhi [velocity vx vy vz] block ...

   Bounds are in box units; INF and -INF give unbounded faces.
------------------------------------------------------------------------- */

#include "fix_wall_bounceback.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double BIG = 1.0e20;

/* ---------------------------------------------------------------------- */

static double bound(const char *str, LAMMPS *lmp)
{
  if (strcmp(str, "INF") == 0) return BIG;
  if (strcmp(str, "-INF") == 0) return -BIG;
  return utils::numeric(FLERR, str, false, lmp);
}

/* ---------------------------------------------------------------------- */

FixWallBounceback::FixWallBounceback(LAMMPS *_lmp, int narg, char **arg) :
    Fix(_lmp, narg, arg), nblock(0), slipbit(0), rround(0.0), nbounce(0.0), nfallback(0.0)
{
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;
  dynamic_group_allow = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "block") == 0) {
      if (iarg + 7 > narg) utils::missing_cmd_args(FLERR, "fix wall/bounceback block", error);
      double b[9];
      for (int k = 0; k < 3; k++) {
        b[k] = bound(arg[iarg + 1 + 2 * k], lmp);
        b[3 + k] = bound(arg[iarg + 2 + 2 * k], lmp);
        b[6 + k] = 0.0;
        if (b[k] >= b[3 + k])
          error->all(FLERR, "Fix wall/bounceback block {} has lo >= hi", nblock + 1);
      }
      block.insert(block.end(), b, b + 9);
      nblock++;
      iarg += 7;
    } else if (strcmp(arg[iarg], "slip") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix wall/bounceback slip", error);
      int jgroup = group->find(arg[iarg + 1]);
      if (jgroup < 0) error->all(FLERR, "Fix wall/bounceback slip group {} does not exist", arg[iarg + 1]);
      slipbit = group->bitmask[jgroup];
      iarg += 2;
    } else if (strcmp(arg[iarg], "round") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix wall/bounceback round", error);
      rround = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (rround < 0.0) error->all(FLERR, "Fix wall/bounceback round must be >= 0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "velocity") == 0) {
      if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix wall/bounceback velocity", error);
      if (nblock == 0) error->all(FLERR, "Fix wall/bounceback velocity must follow a block");
      for (int k = 0; k < 3; k++)
        block[9 * (nblock - 1) + 6 + k] = utils::numeric(FLERR, arg[iarg + 1 + k], false, lmp);
      iarg += 4;
    } else
      error->all(FLERR, "Unknown fix wall/bounceback keyword: {}", arg[iarg]);
  }

  if (nblock == 0) error->all(FLERR, "Fix wall/bounceback needs at least one block");
  if (rround > 0.0 && slipbit == 0) error->all(FLERR, "Fix wall/bounceback round needs slip");
  for (int b = 0; b < nblock; b++)
    for (int k = 0; k < 3; k++)
      if (block[9 * b + 3 + k] - block[9 * b + k] < 2.0 * rround)
        error->all(FLERR, "Fix wall/bounceback block {} is thinner than twice the round radius",
                   b + 1);
}

/* ---------------------------------------------------------------------- */

int FixWallBounceback::setmask()
{
  return POST_INTEGRATE;
}

/* ---------------------------------------------------------------------- */

void FixWallBounceback::init()
{
  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Fix wall/bounceback does not support run_style respa");

  // a block cut by a periodic boundary is only seen in the image inside the box;
  // an open (INF) face is fine, it also catches atoms that have just left the box

  if (comm->me == 0) {
    for (int b = 0; b < nblock; b++)
      for (int k = 0; k < 3; k++) {
        double lo = block[9 * b + k], hi = block[9 * b + 3 + k];
        if (!domain->periodicity[k]) continue;
        if ((lo != -BIG && lo < domain->boxlo[k]) || (hi != BIG && hi > domain->boxhi[k]))
          error->warning(FLERR,
                         "Fix wall/bounceback block {} crosses a periodic boundary; only "
                         "its part inside the box is a wall",
                         b + 1);
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallBounceback::setup(int /*vflag*/)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  bigint n = 0, nall;
  for (int i = 0; i < nlocal; i++)
    if ((mask[i] & groupbit) && inwall(x[i], mask[i]) >= 0) n++;
  MPI_Allreduce(&n, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (nall && comm->me == 0)
    error->warning(FLERR, "Fix wall/bounceback: {} atoms of group {} start inside a wall", nall,
                   group->names[igroup]);
}

/* ---------------------------------------------------------------------- */

void FixWallBounceback::post_integrate()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double xold[3], vold[3], n[3], depth;
  int face;

  nbounce = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    const bool rounded = (mask[i] & slipbit) && rround > 0.0;
    int b = rounded ? inside_round(x[i], n, depth) : inside(x[i]);
    if (b < 0) continue;

    for (int k = 0; k < 3; k++) {
      vold[k] = v[i][k];
      xold[k] = x[i][k] - dt * v[i][k];
    }
    const double *u = &block[9 * b + 6];
    double s = rounded ? 0.0 : entry(b, xold, x[i], face);
    if (rounded) {
      // mirror at the rounded surface: reverse the normal velocity only
      // (n = 0 if the atom is deeper than R; then the old position is restored)
      double vn = 0.0;
      for (int k = 0; k < 3; k++) vn += (v[i][k] - u[k]) * n[k];
      for (int k = 0; k < 3; k++) {
        x[i][k] += 2.0 * depth * n[k];
        if (vn < 0.0) v[i][k] -= 2.0 * vn * n[k];
      }
    } else if ((mask[i] & slipbit) && face >= 0) {
      // mirror at the face: reverse the normal velocity only
      int k = face % 3;
      double plane = block[9 * b + face];
      x[i][k] = 2.0 * plane - x[i][k];
      v[i][k] = 2.0 * u[k] - v[i][k];
    } else {
      for (int k = 0; k < 3; k++) {
        v[i][k] = 2.0 * u[k] - vold[k];
        x[i][k] = xold[k] + s * (x[i][k] - xold[k]) + (1.0 - s) * dt * v[i][k];
      }
    }
    if (inwall(x[i], mask[i]) >= 0) {
      for (int k = 0; k < 3; k++) {
        x[i][k] = xold[k];
        v[i][k] = 2.0 * u[k] - vold[k];
      }
      nfallback += 1.0;
    }
    nbounce += 1.0;
  }
}

/* ----------------------------------------------------------------------
   0 = number of reflections in the last step
   1 = number of reflections so far that restored the old position
------------------------------------------------------------------------- */

double FixWallBounceback::compute_vector(int n)
{
  double one = (n == 0) ? nbounce : nfallback;
  double all;
  MPI_Allreduce(&one, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/* ----------------------------------------------------------------------
   index of the first block that contains point p, -1 if none
------------------------------------------------------------------------- */

int FixWallBounceback::inside(const double *p) const
{
  for (int b = 0; b < nblock; b++) {
    const double *lo = &block[9 * b];
    const double *hi = &block[9 * b + 3];
    if (p[0] > lo[0] && p[0] < hi[0] && p[1] > lo[1] && p[1] < hi[1] && p[2] > lo[2] &&
        p[2] < hi[2])
      return b;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   index of the first block that contains point p when the edges of the
   blocks are rounded with radius rround, -1 if none.  For that block, n
   is the outward normal of the rounded surface nearest to p and depth
   the distance of p below it; n = 0 if p is deeper than rround.
------------------------------------------------------------------------- */

int FixWallBounceback::inside_round(const double *p, double *n, double &depth) const
{
  for (int b = 0; b < nblock; b++) {
    const double *lo = &block[9 * b];
    const double *hi = &block[9 * b + 3];
    double d[3], r2 = 0.0;
    for (int k = 0; k < 3; k++) {
      d[k] = p[k] - MIN(MAX(p[k], lo[k] + rround), hi[k] - rround);
      r2 += d[k] * d[k];
    }
    if (r2 >= rround * rround) continue;
    const double r = sqrt(r2);
    for (int k = 0; k < 3; k++) n[k] = (r > 0.0) ? d[k] / r : 0.0;
    depth = rround - r;
    return b;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   the wall test that applies to an atom with this mask
------------------------------------------------------------------------- */

int FixWallBounceback::inwall(const double *p, int m) const
{
  double n[3], depth;
  if ((m & slipbit) && rround > 0.0) return inside_round(p, n, depth);
  return inside(p);
}

/* ----------------------------------------------------------------------
   fraction of the step xold -> xnew after which the path enters block b,
   and the face it crosses there (0-2 = lower x, y, z face, 3-5 = upper
   faces, -1 if xold is already inside the block)
------------------------------------------------------------------------- */

double FixWallBounceback::entry(int b, const double *xold, const double *xnew, int &face) const
{
  const double *lo = &block[9 * b];
  const double *hi = &block[9 * b + 3];
  double s = 0.0, t;
  face = -1;
  for (int k = 0; k < 3; k++) {
    if (xold[k] <= lo[k]) {
      t = (lo[k] - xold[k]) / (xnew[k] - xold[k]);
      if (face < 0 || t > s) {
        s = t;
        face = k;
      }
    } else if (xold[k] >= hi[k]) {
      t = (hi[k] - xold[k]) / (xnew[k] - xold[k]);
      if (face < 0 || t > s) {
        s = t;
        face = 3 + k;
      }
    }
  }
  return MIN(MAX(s, 0.0), 1.0);
}
