/* ----------------------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.

   Bending energy between two neighboring triangles of the membrane mesh
   (Fedosov, Caswell and Karniadakis, Biophys. J. 98, 2215 (2010)):

     E = kb [1 - cos(theta - theta0)]

   theta is the angle between the normals of triangles 1-2-3 and 4-3-2,
   which share the edge 2-3.  Both triangles are listed counterclockwise
   when seen from outside.  The forces follow dihedral style bend of the
   in-house code lmp_rbc.
------------------------------------------------------------------------- */

#include "dihedral_bend.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::RAD2DEG;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

DihedralBend::DihedralBend(LAMMPS *_lmp) :
    Dihedral(_lmp), k(nullptr), theta0(nullptr), cos0(nullptr), sin0(nullptr)
{
  // the per-atom data compute() reads and writes (with KOKKOS, these are
  // copied between device and host every step)
  datamask_read = X_MASK | F_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

DihedralBend::~DihedralBend()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(theta0);
    memory->destroy(cos0);
    memory->destroy(sin0);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBend::compute(int eflag, int vflag)
{
  int i1, i2, i3, i4, n, type;
  double d21x, d21y, d21z, d31x, d31y, d31z, d32x, d32y, d32z;
  double d34x, d34y, d34z, d24x, d24y, d24z, d14x, d14y, d14z;
  double n1x, n1y, n1z, n2x, n2y, n2z, n1, n2, nn;
  double costheta, sintheta, mx, alfa, a11, a12, a22;
  double f1[3], f2[3], f3[3], f4[3], edihedral;

  edihedral = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    d21x = x[i2][0] - x[i1][0];
    d21y = x[i2][1] - x[i1][1];
    d21z = x[i2][2] - x[i1][2];
    d31x = x[i3][0] - x[i1][0];
    d31y = x[i3][1] - x[i1][1];
    d31z = x[i3][2] - x[i1][2];
    d32x = x[i3][0] - x[i2][0];
    d32y = x[i3][1] - x[i2][1];
    d32z = x[i3][2] - x[i2][2];
    d34x = x[i3][0] - x[i4][0];
    d34y = x[i3][1] - x[i4][1];
    d34z = x[i3][2] - x[i4][2];
    d24x = x[i2][0] - x[i4][0];
    d24y = x[i2][1] - x[i4][1];
    d24z = x[i2][2] - x[i4][2];
    d14x = x[i1][0] - x[i4][0];
    d14y = x[i1][1] - x[i4][1];
    d14z = x[i1][2] - x[i4][2];

    // normals of the two triangles

    n1x = d21y * d31z - d31y * d21z;
    n1y = d31x * d21z - d21x * d31z;
    n1z = d21x * d31y - d31x * d21y;
    n2x = d34y * d24z - d24y * d34z;
    n2y = d24x * d34z - d34x * d24z;
    n2z = d34x * d24y - d24x * d34y;
    n1 = n1x * n1x + n1y * n1y + n1z * n1z;
    n2 = n2x * n2x + n2y * n2y + n2z * n2z;
    nn = sqrt(n1 * n2);

    // theta > 0 when the surface is convex at the shared edge

    costheta = (n1x * n2x + n1y * n2y + n1z * n2z) / nn;
    if (costheta > 1.0) costheta = 1.0;
    if (costheta < -1.0) costheta = -1.0;
    sintheta = sqrt(1.0 - costheta * costheta);
    if (sintheta < SMALL) sintheta = SMALL;
    mx = (n1x - n2x) * d14x + (n1y - n2y) * d14y + (n1z - n2z) * d14z;
    if (mx < 0.0) sintheta = -sintheta;

    alfa = k[type] * (cos0[type] - costheta * sin0[type] / sintheta);
    a11 = -alfa * costheta / n1;
    a12 = alfa / nn;
    a22 = -alfa * costheta / n2;

    f1[0] = a11 * (n1y * d32z - n1z * d32y) + a12 * (n2y * d32z - n2z * d32y);
    f1[1] = a11 * (n1z * d32x - n1x * d32z) + a12 * (n2z * d32x - n2x * d32z);
    f1[2] = a11 * (n1x * d32y - n1y * d32x) + a12 * (n2x * d32y - n2y * d32x);
    f2[0] = a11 * (n1z * d31y - n1y * d31z) + a22 * (n2y * d34z - n2z * d34y) +
        a12 * (n2z * d31y - n2y * d31z + n1y * d34z - n1z * d34y);
    f2[1] = a11 * (n1x * d31z - n1z * d31x) + a22 * (n2z * d34x - n2x * d34z) +
        a12 * (n2x * d31z - n2z * d31x + n1z * d34x - n1x * d34z);
    f2[2] = a11 * (n1y * d31x - n1x * d31y) + a22 * (n2x * d34y - n2y * d34x) +
        a12 * (n2y * d31x - n2x * d31y + n1x * d34y - n1y * d34x);
    f3[0] = a11 * (n1y * d21z - n1z * d21y) + a22 * (n2z * d24y - n2y * d24z) +
        a12 * (n2y * d21z - n2z * d21y + n1z * d24y - n1y * d24z);
    f3[1] = a11 * (n1z * d21x - n1x * d21z) + a22 * (n2x * d24z - n2z * d24x) +
        a12 * (n2z * d21x - n2x * d21z + n1x * d24z - n1z * d24x);
    f3[2] = a11 * (n1x * d21y - n1y * d21x) + a22 * (n2y * d24x - n2x * d24y) +
        a12 * (n2x * d21y - n2y * d21x + n1y * d24x - n1x * d24y);
    f4[0] = a22 * (n2z * d32y - n2y * d32z) + a12 * (n1z * d32y - n1y * d32z);
    f4[1] = a22 * (n2x * d32z - n2z * d32x) + a12 * (n1x * d32z - n1z * d32x);
    f4[2] = a22 * (n2y * d32x - n2x * d32y) + a12 * (n1y * d32x - n1x * d32y);

    if (eflag) edihedral = k[type] * (1.0 - (costheta * cos0[type] + sintheta * sin0[type]));

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (evflag)
      ev_tally(i1, i2, i3, i4, nlocal, newton_bond, edihedral, f1, f3, f4, -d21x, -d21y, -d21z,
               d32x, d32y, d32z, -d34x, -d34y, -d34z);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBend::allocate()
{
  allocated = 1;
  const int np1 = atom->ndihedraltypes + 1;

  memory->create(k, np1, "dihedral:k");
  memory->create(theta0, np1, "dihedral:theta0");
  memory->create(cos0, np1, "dihedral:cos0");
  memory->create(sin0, np1, "dihedral:sin0");
  memory->create(setflag, np1, "dihedral:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   dihedral_coeff N kb theta0 (degrees)
------------------------------------------------------------------------- */

void DihedralBend::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Dihedral bend needs 2 coefficients: kb theta0");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->ndihedraltypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double theta0_one = utils::numeric(FLERR, arg[2], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    theta0[i] = theta0_one * DEG2RAD;
    cos0[i] = cos(theta0[i]);
    sin0[i] = sin(theta0[i]);
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for dihedral coefficients");
}

/* ---------------------------------------------------------------------- */

void DihedralBend::write_restart(FILE *fp)
{
  fwrite(&k[1], sizeof(double), atom->ndihedraltypes, fp);
  fwrite(&theta0[1], sizeof(double), atom->ndihedraltypes, fp);
}

/* ---------------------------------------------------------------------- */

void DihedralBend::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->ndihedraltypes, fp, nullptr, error);
    utils::sfread(FLERR, &theta0[1], sizeof(double), atom->ndihedraltypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta0[1], atom->ndihedraltypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    cos0[i] = cos(theta0[i]);
    sin0[i] = sin(theta0[i]);
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBend::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    fprintf(fp, "%d %g %g\n", i, k[i], theta0[i] * RAD2DEG);
}
