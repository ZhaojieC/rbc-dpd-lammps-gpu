/* ----------------------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.

   If you use this code, please cite Z. Chai, M. de Person, P. A. Buffet,
   M. Franco and G. E. Karniadakis, Biophys. J. (2026),
   https://doi.org/10.1016/j.bpj.2026.10.021.  LAMMPS prints this
   reminder at the end of a run and writes the BibTeX entries to the
   log file (see rbc_cite.h).

   compute ID all rbc/shape

   Global array with one row per molecule ID (row N = molecule N):
     1 = membrane area, 2 = enclosed volume,
     3 = reduced volume 6 sqrt(pi) V / A^(3/2) (1 for a sphere).
   Rows of molecules that are not cells are zero.  The values come from
   angle style area/volume, which must be in use.
------------------------------------------------------------------------- */

#include "compute_rbc_shape.h"

#include "angle_area_volume.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "rbc_cite.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

ComputeRBCShape::ComputeRBCShape(LAMMPS *_lmp, int narg, char **arg) :
    Compute(_lmp, narg, arg), nrows(0), data(nullptr), av(nullptr)
{
  rbc_cite(lmp, false);    // the articles to cite, printed after the run

  if (narg != 3) error->all(FLERR, "Illegal compute rbc/shape command");
  if (!atom->molecule_flag) error->all(FLERR, "Compute rbc/shape requires molecule IDs");

  tagint mymax = 0, allmax;
  for (int i = 0; i < atom->nlocal; i++) mymax = MAX(mymax, atom->molecule[i]);
  MPI_Allreduce(&mymax, &allmax, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  if (allmax < 1) error->all(FLERR, "Compute rbc/shape found no molecules");
  nrows = allmax;

  array_flag = 1;
  size_array_rows = nrows;
  size_array_cols = 3;
  extarray = 0;

  memory->create(data, nrows, 3, "rbc/shape:data");
  array = data;
}

/* ---------------------------------------------------------------------- */

ComputeRBCShape::~ComputeRBCShape()
{
  memory->destroy(data);
}

/* ---------------------------------------------------------------------- */

void ComputeRBCShape::init()
{
  av = dynamic_cast<AngleAreaVolume *>(force->angle);
  if (!av) error->all(FLERR, "Compute rbc/shape requires angle style area/volume");
}

/* ---------------------------------------------------------------------- */

void ComputeRBCShape::compute_array()
{
  invoked_array = update->ntimestep;

  double *area, *volume;
  int *celltype;
  int maxmol = av->cell_totals(area, volume, celltype);

  for (int m = 1; m <= nrows; m++) {
    double *row = data[m - 1];
    row[0] = row[1] = row[2] = 0.0;
    if (m > maxmol || celltype[m] <= 0) continue;
    row[0] = area[m];
    row[1] = volume[m];
    if (area[m] > 0.0) row[2] = 6.0 * sqrt(MY_PI) * volume[m] / pow(area[m], 1.5);
  }
}
