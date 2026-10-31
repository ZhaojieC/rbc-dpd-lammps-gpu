/* ----------------------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.

   WLC-POW spring of the RBC spectrin network (Fedosov, Caswell and
   Karniadakis, Biophys. J. 98, 2215 (2010)).  Rewritten for current
   LAMMPS from the rbc branch of bond style break/ligand in the in-house
   code lmp_rbc (contributing authors Alireza Yazdani and Dmitry Fedosov).

   Each bond keeps its own equilibrium length l0, taken from the reference
   shape (see rbc_reference.h).  For every bond

     lmax = ratio * l0
     F(r) = -(1/lambda) [1/(4(1-x)^2) - 1/4 + x] + kp / r^m,   x = r / lmax

   where lambda and kp make F(l0) = 0 and give the network the linear
   shear modulus mu0 (Fedosov et al. 2010).
------------------------------------------------------------------------- */

#include "bond_wlc_pow.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "rbc_reference.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondWLCPOW::BondWLCPOW(LAMMPS *_lmp) :
    Bond(_lmp), ratio(nullptr), mu0(nullptr), expo(nullptr), ref(nullptr)
{
  // the per-atom data compute() reads and writes (with KOKKOS, these are
  // copied between device and host every step)
  datamask_read = X_MASK | F_MASK | TAG_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

BondWLCPOW::~BondWLCPOW()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(ratio);
    memory->destroy(mu0);
    memory->destroy(expo);
  }
}

/* ---------------------------------------------------------------------- */

void BondWLCPOW::compute(int eflag, int vflag)
{
  double lmax, lambda, kp;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (int n = 0; n < nbondlist; n++) {
    int i1 = bondlist[n][0];
    int i2 = bondlist[n][1];
    int type = bondlist[n][2];

    double delx = x[i1][0] - x[i2][0];
    double dely = x[i1][1] - x[i2][1];
    double delz = x[i1][2] - x[i2][2];
    double r = sqrt(delx * delx + dely * dely + delz * delz);

    spring(type, reference_length(i1, i2), lmax, lambda, kp);
    double s = r / lmax;
    if (s >= 1.0)
      error->one(FLERR, "Bond {}-{} reached its maximum length at step {} (r/lmax = {:.4})",
                 atom->tag[i1], atom->tag[i2], update->ntimestep, s);

    double m = expo[type];
    double fbond = -(0.25 / ((1.0 - s) * (1.0 - s)) - 0.25 + s) / lambda / r + kp / pow(r, m + 1.0);

    double ebond = 0.0;
    if (eflag) {
      ebond = 0.25 * lmax * (3.0 * s * s - 2.0 * s * s * s) / (1.0 - s) / lambda;
      if (m == 1.0)
        ebond -= kp * log(r);
      else
        ebond += kp / (m - 1.0) / pow(r, m - 1.0);
    }

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx * fbond;
      f[i1][1] += dely * fbond;
      f[i1][2] += delz * fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx * fbond;
      f[i2][1] -= dely * fbond;
      f[i2][2] -= delz * fbond;
    }

    if (evflag) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
  }
}

/* ----------------------------------------------------------------------
   length of the bond i-j in the reference shape
------------------------------------------------------------------------- */

double BondWLCPOW::reference_length(int i, int j)
{
  double **x0 = ref->astore;
  double dx = x0[i][0] - x0[j][0];
  double dy = x0[i][1] - x0[j][1];
  double dz = x0[i][2] - x0[j][2];
  domain->minimum_image(FLERR, dx, dy, dz);
  return sqrt(dx * dx + dy * dy + dz * dz);
}

/* ----------------------------------------------------------------------
   lmax, 1/lambda and kp of a bond with equilibrium length l0
   (same algebra as lmp_rbc; the thermal energy of the original
   coefficients cancels and is not needed)
------------------------------------------------------------------------- */

void BondWLCPOW::spring(int type, double l0, double &lmax, double &lambda, double &kp)
{
  const double m = expo[type];
  const double rr = 1.0 / ratio[type];

  lmax = l0 * ratio[type];
  kp = pow(l0, m) * (0.25 / ((1.0 - rr) * (1.0 - rr)) - 0.25 + rr);
  double mu = 0.25 * sqrt(3.0) *
      ((-0.25 / pow(1.0 - rr, 2.0) + 0.25 + 0.5 * rr / pow(1.0 - rr, 3.0)) / lmax / rr +
       kp * (m + 1.0) / pow(l0, m + 1.0));
  lambda = mu / mu0[type];
  kp *= mu0[type] / mu;
}

/* ---------------------------------------------------------------------- */

void BondWLCPOW::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(ratio, np1, "bond:ratio");
  memory->create(mu0, np1, "bond:mu0");
  memory->create(expo, np1, "bond:expo");
  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   bond_coeff N lmax/l0 mu0 m
------------------------------------------------------------------------- */

void BondWLCPOW::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR, "Bond wlc/pow needs 3 coefficients: lmax/l0 mu0 m");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double ratio_one = utils::numeric(FLERR, arg[1], false, lmp);
  double mu0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double expo_one = utils::numeric(FLERR, arg[3], false, lmp);
  if (ratio_one <= 1.0) error->all(FLERR, "Bond wlc/pow lmax/l0 must be larger than 1");
  if (mu0_one <= 0.0) error->all(FLERR, "Bond wlc/pow mu0 must be positive");
  if (expo_one <= 0.0) error->all(FLERR, "Bond wlc/pow exponent m must be positive");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    ratio[i] = ratio_one;
    mu0[i] = mu0_one;
    expo[i] = expo_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ---------------------------------------------------------------------- */

void BondWLCPOW::init_style()
{
  ref = rbc_reference(lmp);
}

/* ---------------------------------------------------------------------- */

void BondWLCPOW::write_restart(FILE *fp)
{
  fwrite(&ratio[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&mu0[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&expo[1], sizeof(double), atom->nbondtypes, fp);
}

/* ---------------------------------------------------------------------- */

void BondWLCPOW::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &ratio[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &mu0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &expo[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&ratio[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&mu0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&expo[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

void BondWLCPOW::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp, "%d %g %g %g\n", i, ratio[i], mu0[i], expo[i]);
}

/* ---------------------------------------------------------------------- */

double BondWLCPOW::single(int type, double rsq, int i, int j, double &fforce)
{
  double lmax, lambda, kp;
  spring(type, reference_length(i, j), lmax, lambda, kp);

  double r = sqrt(rsq);
  double s = r / lmax;
  double m = expo[type];
  fforce = -(0.25 / ((1.0 - s) * (1.0 - s)) - 0.25 + s) / lambda / r + kp / pow(r, m + 1.0);

  double eng = 0.25 * lmax * (3.0 * s * s - 2.0 * s * s * s) / (1.0 - s) / lambda;
  if (m == 1.0)
    eng -= kp * log(r);
  else
    eng += kp / (m - 1.0) / pow(r, m - 1.0);
  return eng;
}
