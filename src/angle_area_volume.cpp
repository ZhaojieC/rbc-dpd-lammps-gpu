/* ----------------------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.

   If you use this code, please cite Z. Chai, M. de Person, P. A. Buffet,
   M. Franco and G. E. Karniadakis, Biophys. J. (2026),
   https://doi.org/10.1016/j.bpj.2026.10.021.  LAMMPS prints this
   reminder at the end of a run and writes the BibTeX entries to the
   log file (see rbc_cite.h).

   Area and volume constraints of a triangulated cell membrane (Fedosov,
   Caswell and Karniadakis, Biophys. J. 98, 2215 (2010)).  Each angle is
   one triangle of the mesh, listed counterclockwise when seen from
   outside.  Each cell is one molecule.  The energy is

     E = sum_cells [ ka (A - A0)^2 / (2 A0) + kv (V - V0)^2 / (2 V0) ]
       + sum_triangles kd (a - a0)^2 / (2 a0)

   with A and V the total area and enclosed volume of the cell and a0
   the area of the triangle in the reference shape (rbc_reference.h).
   The forces follow angle style area/volume of the in-house code
   lmp_rbc.  Here the volume is computed relative to a reference point of
   each cell, so a cell may cross a periodic boundary.
------------------------------------------------------------------------- */

#include "angle_area_volume.h"

#include "accelerator_kokkos.h"
#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "rbc_cite.h"
#include "rbc_reference.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AngleAreaVolume::AngleAreaVolume(LAMMPS *_lmp) :
    Angle(_lmp), ka(nullptr), area0(nullptr), kv(nullptr), volume0(nullptr), kd(nullptr),
    ref(nullptr), maxmol(-1), celltype(nullptr), nmolatoms(nullptr), com(nullptr), area(nullptr),
    volume(nullptr), buf1(nullptr), buf2(nullptr)
{
  rbc_cite(lmp, true);    // the articles to cite, printed after the run

  // the volume term is not a sum of per-triangle virials
  centroidstressflag = CENTROID_NOTAVAIL;

  // the per-atom data compute() reads and writes (with KOKKOS, these are
  // copied between device and host every step)
  datamask_read = X_MASK | F_MASK | IMAGE_MASK | MOLECULE_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

AngleAreaVolume::~AngleAreaVolume()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(ka);
    memory->destroy(area0);
    memory->destroy(kv);
    memory->destroy(volume0);
    memory->destroy(kd);
  }
  memory->destroy(celltype);
  memory->destroy(nmolatoms);
  memory->destroy(com);
  memory->destroy(area);
  memory->destroy(volume);
  memory->destroy(buf1);
  memory->destroy(buf2);
}

/* ---------------------------------------------------------------------- */

void AngleAreaVolume::compute(int eflag, int vflag)
{
  int i1, i2, i3, n, m, type;
  double d21x, d21y, d21z, d31x, d31y, d31z, d32x, d32y, d32z;
  double nx, ny, nz, nn, mx, my, mz, r1x, r1y, r1z;
  double a0, coefl, coefa, coefca, coefv, eangle;
  double s1[3], s2[3], s3[3], v1[3], v2[3], v3[3];

  ev_init(eflag, vflag);

  reference_points();
  totals();

  double **x = atom->x;
  double **f = atom->f;
  tagint *molecule = atom->molecule;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
    m = molecule[i2];

    d21x = x[i2][0] - x[i1][0];
    d21y = x[i2][1] - x[i1][1];
    d21z = x[i2][2] - x[i1][2];
    d31x = x[i3][0] - x[i1][0];
    d31y = x[i3][1] - x[i1][1];
    d31z = x[i3][2] - x[i1][2];
    d32x = x[i3][0] - x[i2][0];
    d32y = x[i3][1] - x[i2][1];
    d32z = x[i3][2] - x[i2][2];

    // normal vector, its length is twice the triangle area

    nx = d21y * d31z - d31y * d21z;
    ny = d31x * d21z - d21x * d31z;
    nz = d21x * d31y - d31x * d21y;
    nn = sqrt(nx * nx + ny * ny + nz * nz);

    // sum of the three vertices relative to the reference point of the cell

    r1x = x[i1][0] - com[3 * m];
    r1y = x[i1][1] - com[3 * m + 1];
    r1z = x[i1][2] - com[3 * m + 2];
    domain->minimum_image_big(FLERR, r1x, r1y, r1z);
    mx = 3.0 * r1x + d21x + d31x;
    my = 3.0 * r1y + d21y + d31y;
    mz = 3.0 * r1z + d21z + d31z;

    a0 = triangle_area0(i1, i2, i3);
    coefl = 0.5 * kd[type] * (a0 - 0.5 * nn) / a0 / nn;
    coefa = 0.5 * ka[type] * (area0[type] - area[m]) / area0[type] / nn;
    coefca = coefl + coefa;
    coefv = kv[type] * (volume0[type] - volume[m]) / volume0[type] / 18.0;

    // area forces, they add up to zero on each triangle

    s1[0] = coefca * (ny * d32z - nz * d32y);
    s1[1] = coefca * (nz * d32x - nx * d32z);
    s1[2] = coefca * (nx * d32y - ny * d32x);
    s2[0] = coefca * (nz * d31y - ny * d31z);
    s2[1] = coefca * (nx * d31z - nz * d31x);
    s2[2] = coefca * (ny * d31x - nx * d31y);
    s3[0] = coefca * (ny * d21z - nz * d21y);
    s3[1] = coefca * (nz * d21x - nx * d21z);
    s3[2] = coefca * (nx * d21y - ny * d21x);

    // volume forces

    v1[0] = coefv * (nx + d32z * my - d32y * mz);
    v1[1] = coefv * (ny - d32z * mx + d32x * mz);
    v1[2] = coefv * (nz + d32y * mx - d32x * my);
    v2[0] = coefv * (nx - d31z * my + d31y * mz);
    v2[1] = coefv * (ny + d31z * mx - d31x * mz);
    v2[2] = coefv * (nz - d31y * mx + d31x * my);
    v3[0] = coefv * (nx + d21z * my - d21y * mz);
    v3[1] = coefv * (ny - d21z * mx + d21x * mz);
    v3[2] = coefv * (nz + d21y * mx - d21x * my);

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += s1[0] + v1[0];
      f[i1][1] += s1[1] + v1[1];
      f[i1][2] += s1[2] + v1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += s2[0] + v2[0];
      f[i2][1] += s2[1] + v2[1];
      f[i2][2] += s2[2] + v2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += s3[0] + v3[0];
      f[i3][1] += s3[1] + v3[1];
      f[i3][2] += s3[2] + v3[2];
    }

    // the per-triangle tally covers the local area energy and the virial
    // of the area forces; the cell terms are added below

    if (evflag) {
      eangle = 0.0;
      if (eflag) eangle = 0.5 * kd[type] * (a0 - 0.5 * nn) * (a0 - 0.5 * nn) / a0;
      ev_tally(i1, i2, i3, nlocal, newton_bond, eangle, s1, s3, -d21x, -d21y, -d21z, d32x, d32y,
               d32z);
    }
  }

  if (!evflag) return;

  // energy of the total-area and volume terms of each cell, and the virial
  // of the volume forces: for a closed surface sum_i x_i dV/dx_i = V I,
  // so that virial is isotropic.  buf1 = energy, buf2 = virial per cell.

  for (m = 0; m <= maxmol; m++) {
    buf1[m] = buf2[m] = 0.0;
    type = celltype[m];
    if (type <= 0) continue;
    double da = area[m] - area0[type];
    double dv = volume[m] - volume0[type];
    buf1[m] = 0.5 * ka[type] * da * da / area0[type] + 0.5 * kv[type] * dv * dv / volume0[type];
    buf2[m] = -kv[type] * dv * volume[m] / volume0[type];
    if (comm->me == 0) {
      if (eflag_global) energy += buf1[m];
      if (vflag_global) {
        virial[0] += buf2[m];
        virial[1] += buf2[m];
        virial[2] += buf2[m];
      }
    }
  }

  if (eflag_atom || vflag_atom) {
    for (int i = 0; i < nlocal; i++) {
      m = molecule[i];
      if (m <= 0 || m > maxmol || celltype[m] <= 0) continue;
      if (eflag_atom) eatom[i] += buf1[m] / nmolatoms[m];
      if (vflag_atom) {
        vatom[i][0] += buf2[m] / nmolatoms[m];
        vatom[i][1] += buf2[m] / nmolatoms[m];
        vatom[i][2] += buf2[m] / nmolatoms[m];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reference point of each cell: its center of mass (equal weights)
------------------------------------------------------------------------- */

void AngleAreaVolume::reference_points()
{
  double **x = atom->x;
  imageint *image = atom->image;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nbuf = 4 * (maxmol + 1);
  double xu[3];

  for (int k = 0; k < nbuf; k++) buf1[k] = 0.0;

  for (int i = 0; i < nlocal; i++) {
    tagint m = molecule[i];
    if (m <= 0 || m > maxmol || celltype[m] <= 0) continue;
    domain->unmap(x[i], image[i], xu);
    buf1[4 * m] += xu[0];
    buf1[4 * m + 1] += xu[1];
    buf1[4 * m + 2] += xu[2];
    buf1[4 * m + 3] += 1.0;
  }

  MPI_Allreduce(buf1, buf2, nbuf, MPI_DOUBLE, MPI_SUM, world);

  for (int m = 0; m <= maxmol; m++) {
    if (celltype[m] <= 0) continue;
    com[3 * m] = buf2[4 * m] / buf2[4 * m + 3];
    com[3 * m + 1] = buf2[4 * m + 1] / buf2[4 * m + 3];
    com[3 * m + 2] = buf2[4 * m + 2] / buf2[4 * m + 3];
  }

  // the volume is summed relative to the center, so every vertex must lie
  // within half a box length of it

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    tagint m = molecule[i];
    if (m <= 0 || m > maxmol || celltype[m] <= 0) continue;
    domain->unmap(x[i], image[i], xu);
    for (int k = 0; k < 3; k++)
      if (domain->periodicity[k] && fabs(xu[k] - com[3 * m + k]) > 0.5 * domain->prd[k]) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
  if (flagall)
    error->all(FLERR,
               "Angle style area/volume: a cell is larger than half the periodic box, or its "
               "image flags are inconsistent");
}

/* ----------------------------------------------------------------------
   total area and enclosed volume of each cell
   a triangle is counted on the processor that owns its second atom
------------------------------------------------------------------------- */

void AngleAreaVolume::totals()
{
  double **x = atom->x;
  tagint *molecule = atom->molecule;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int nbuf = 2 * (maxmol + 1);

  for (int k = 0; k < nbuf; k++) buf1[k] = 0.0;

  for (int n = 0; n < nanglelist; n++) {
    int i1 = anglelist[n][0];
    int i2 = anglelist[n][1];
    int i3 = anglelist[n][2];
    if (i2 >= nlocal) continue;
    int m = molecule[i2];

    double d21x = x[i2][0] - x[i1][0];
    double d21y = x[i2][1] - x[i1][1];
    double d21z = x[i2][2] - x[i1][2];
    double d31x = x[i3][0] - x[i1][0];
    double d31y = x[i3][1] - x[i1][1];
    double d31z = x[i3][2] - x[i1][2];

    double nx = d21y * d31z - d31y * d21z;
    double ny = d31x * d21z - d21x * d31z;
    double nz = d21x * d31y - d31x * d21y;

    double r1x = x[i1][0] - com[3 * m];
    double r1y = x[i1][1] - com[3 * m + 1];
    double r1z = x[i1][2] - com[3 * m + 2];
    domain->minimum_image_big(FLERR, r1x, r1y, r1z);
    double mx = 3.0 * r1x + d21x + d31x;
    double my = 3.0 * r1y + d21y + d31y;
    double mz = 3.0 * r1z + d21z + d31z;

    buf1[2 * m] += 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
    buf1[2 * m + 1] += (nx * mx + ny * my + nz * mz) / 18.0;
  }

  MPI_Allreduce(buf1, buf2, nbuf, MPI_DOUBLE, MPI_SUM, world);

  for (int m = 0; m <= maxmol; m++) {
    area[m] = buf2[2 * m];
    volume[m] = buf2[2 * m + 1];
    int type = celltype[m];
    if (type > 0 && kv[type] != 0.0 && volume[m] <= 0.0)
      error->all(FLERR,
                 "Cell (molecule {}) has volume {:.6g}; its triangles must be listed "
                 "counterclockwise when seen from outside",
                 m, volume[m]);
  }
}

/* ----------------------------------------------------------------------
   area of triangle i1-i2-i3 in the reference shape
------------------------------------------------------------------------- */

double AngleAreaVolume::triangle_area0(int i1, int i2, int i3)
{
  double **x0 = ref->astore;
  double ax = x0[i2][0] - x0[i1][0];
  double ay = x0[i2][1] - x0[i1][1];
  double az = x0[i2][2] - x0[i1][2];
  double bx = x0[i3][0] - x0[i1][0];
  double by = x0[i3][1] - x0[i1][1];
  double bz = x0[i3][2] - x0[i1][2];
  domain->minimum_image(FLERR, ax, ay, az);
  domain->minimum_image(FLERR, bx, by, bz);
  double nx = ay * bz - by * az;
  double ny = bx * az - ax * bz;
  double nz = ax * by - bx * ay;
  return 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
}

/* ----------------------------------------------------------------------
   find the cells (molecules with angles) and check them
------------------------------------------------------------------------- */

void AngleAreaVolume::setup_cells()
{
  if (!atom->molecule_flag) error->all(FLERR, "Angle style area/volume requires molecule IDs");
  if (lmp->kokkos)
    dynamic_cast<AtomKokkos *>(atom)->sync(Host, TAG_MASK | MOLECULE_MASK | ANGLE_MASK);

  tagint *molecule = atom->molecule;
  tagint *tag = atom->tag;
  int *num_angle = atom->num_angle;
  int **angle_type = atom->angle_type;
  tagint **angle_atom2 = atom->angle_atom2;
  int nlocal = atom->nlocal;

  tagint mymax = 0, allmax;
  for (int i = 0; i < nlocal; i++) mymax = MAX(mymax, molecule[i]);
  MPI_Allreduce(&mymax, &allmax, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  if (allmax > MAXSMALLINT / 4) error->all(FLERR, "Angle style area/volume: molecule IDs too large");

  if (allmax != maxmol) {
    maxmol = allmax;
    memory->destroy(celltype);
    memory->destroy(nmolatoms);
    memory->destroy(com);
    memory->destroy(area);
    memory->destroy(volume);
    memory->destroy(buf1);
    memory->destroy(buf2);
    memory->create(celltype, maxmol + 1, "angle:celltype");
    memory->create(nmolatoms, maxmol + 1, "angle:nmolatoms");
    memory->create(com, 3 * (maxmol + 1), "angle:com");
    memory->create(area, maxmol + 1, "angle:area");
    memory->create(volume, maxmol + 1, "angle:volume");
    memory->create(buf1, 4 * (maxmol + 1), "angle:buf1");
    memory->create(buf2, 4 * (maxmol + 1), "angle:buf2");
  }

  int *tmin, *tmax, *tall;
  double *count, *countall;
  memory->create(tmin, maxmol + 1, "angle:tmin");
  memory->create(tmax, maxmol + 1, "angle:tmax");
  memory->create(tall, maxmol + 1, "angle:tall");
  memory->create(count, 2 * (maxmol + 1), "angle:count");
  memory->create(countall, 2 * (maxmol + 1), "angle:countall");
  for (int m = 0; m <= maxmol; m++) {
    tmin[m] = MAXSMALLINT;
    tmax[m] = 0;
    count[2 * m] = count[2 * m + 1] = 0.0;
  }

  // count[2m] = atoms, count[2m+1] = triangles (each counted once, by its second atom)

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    int m = molecule[i];
    count[2 * m] += 1.0;
    for (int k = 0; k < num_angle[i]; k++) {
      int t = angle_type[i][k];
      if (t <= 0) continue;
      if (m == 0) flag = 1;
      tmin[m] = MIN(tmin[m], t);
      tmax[m] = MAX(tmax[m], t);
      if (angle_atom2[i][k] == tag[i]) count[2 * m + 1] += 1.0;
    }
  }

  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
  if (flagall) error->all(FLERR, "Angle style area/volume: atoms with angles need a molecule ID > 0");

  MPI_Allreduce(tmin, tall, maxmol + 1, MPI_INT, MPI_MIN, world);
  memcpy(tmin, tall, (maxmol + 1) * sizeof(int));
  MPI_Allreduce(tmax, tall, maxmol + 1, MPI_INT, MPI_MAX, world);
  MPI_Allreduce(count, countall, 2 * (maxmol + 1), MPI_DOUBLE, MPI_SUM, world);

  for (int m = 0; m <= maxmol; m++) {
    celltype[m] = tall[m];
    nmolatoms[m] = countall[2 * m];
    if (tall[m] == 0) continue;
    if (tmin[m] != tall[m])
      error->all(FLERR, "Angle style area/volume: all triangles of molecule {} need the same type",
                 m);
    int type = tall[m];
    double nvert = countall[2 * m];
    double ntri = countall[2 * m + 1];
    if (kv[type] != 0.0 && ntri != 2.0 * nvert - 4.0)
      error->all(FLERR,
                 "Angle style area/volume: molecule {} is not a closed surface "
                 "({} triangles for {} vertices)",
                 m, ntri, nvert);
  }

  memory->destroy(tmin);
  memory->destroy(tmax);
  memory->destroy(tall);
  memory->destroy(count);
  memory->destroy(countall);
}

/* ---------------------------------------------------------------------- */

void AngleAreaVolume::allocate()
{
  allocated = 1;
  const int np1 = atom->nangletypes + 1;

  memory->create(ka, np1, "angle:ka");
  memory->create(area0, np1, "angle:area0");
  memory->create(kv, np1, "angle:kv");
  memory->create(volume0, np1, "angle:volume0");
  memory->create(kd, np1, "angle:kd");
  memory->create(setflag, np1, "angle:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   angle_coeff N ka A0 kv V0 kd
------------------------------------------------------------------------- */

void AngleAreaVolume::coeff(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR, "Angle area/volume needs 5 coefficients: ka A0 kv V0 kd");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  double ka_one = utils::numeric(FLERR, arg[1], false, lmp);
  double area0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double kv_one = utils::numeric(FLERR, arg[3], false, lmp);
  double volume0_one = utils::numeric(FLERR, arg[4], false, lmp);
  double kd_one = utils::numeric(FLERR, arg[5], false, lmp);
  if (area0_one <= 0.0 || volume0_one <= 0.0)
    error->all(FLERR, "Angle area/volume A0 and V0 must be positive");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    ka[i] = ka_one;
    area0[i] = area0_one;
    kv[i] = kv_one;
    volume0[i] = volume0_one;
    kd[i] = kd_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

void AngleAreaVolume::init_style()
{
  if (domain->triclinic) error->all(FLERR, "Angle style area/volume does not support triclinic boxes");
  ref = rbc_reference(lmp);
  setup_cells();
}

/* ---------------------------------------------------------------------- */

void AngleAreaVolume::write_restart(FILE *fp)
{
  fwrite(&ka[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&area0[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&kv[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&volume0[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&kd[1], sizeof(double), atom->nangletypes, fp);
}

/* ---------------------------------------------------------------------- */

void AngleAreaVolume::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &ka[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &area0[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &kv[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &volume0[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &kd[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
  }
  MPI_Bcast(&ka[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&area0[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&kv[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&volume0[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&kd[1], atom->nangletypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

void AngleAreaVolume::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp, "%d %g %g %g %g %g\n", i, ka[i], area0[i], kv[i], volume0[i], kd[i]);
}

/* ----------------------------------------------------------------------
   local area energy of one triangle (the cell terms are not per triangle)
------------------------------------------------------------------------- */

double AngleAreaVolume::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;
  double ax = x[i2][0] - x[i1][0];
  double ay = x[i2][1] - x[i1][1];
  double az = x[i2][2] - x[i1][2];
  double bx = x[i3][0] - x[i1][0];
  double by = x[i3][1] - x[i1][1];
  double bz = x[i3][2] - x[i1][2];
  domain->minimum_image(FLERR, ax, ay, az);
  domain->minimum_image(FLERR, bx, by, bz);
  double nx = ay * bz - by * az;
  double ny = bx * az - ax * bz;
  double nz = ax * by - bx * ay;
  double a = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
  double a0 = triangle_area0(i1, i2, i3);
  return 0.5 * kd[type] * (a - a0) * (a - a0) / a0;
}

/* ---------------------------------------------------------------------- */

void *AngleAreaVolume::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "ka") == 0) return (void *) ka;
  if (strcmp(str, "A0") == 0) return (void *) area0;
  if (strcmp(str, "kv") == 0) return (void *) kv;
  if (strcmp(str, "V0") == 0) return (void *) volume0;
  if (strcmp(str, "kd") == 0) return (void *) kd;
  return nullptr;
}

/* ---------------------------------------------------------------------- */

int AngleAreaVolume::cell_totals(double *&a, double *&v, int *&type)
{
  reference_points();
  totals();
  a = area;
  v = volume;
  type = celltype;
  return maxmol;
}

/* ---------------------------------------------------------------------- */

double AngleAreaVolume::memory_usage()
{
  double bytes = Angle::memory_usage();
  bytes += 5.0 * (atom->nangletypes + 1) * sizeof(double);
  bytes += (maxmol + 1) * (sizeof(int) + 14.0 * sizeof(double));
  return bytes;
}
