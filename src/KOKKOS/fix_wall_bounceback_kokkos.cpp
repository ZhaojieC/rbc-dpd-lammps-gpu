/* ----------------------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps-gpu
   Distributed under the GNU General Public License, version 2.

   If you use this code, please cite Z. Chai, M. de Person, P. A. Buffet,
   M. Franco and G. E. Karniadakis, Biophys. J. (2026),
   https://doi.org/10.1016/j.bpj.2026.10.021.  LAMMPS prints this
   reminder at the end of a run and writes the BibTeX entries to the
   log file (see rbc_cite.h).

   KOKKOS version of fix wall/bounceback: the same reflection, done on
   the device.  The arithmetic is in double precision.
------------------------------------------------------------------------- */

#include "fix_wall_bounceback_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallBouncebackKokkos<DeviceType>::FixWallBouncebackKokkos(LAMMPS *_lmp, int narg, char **arg) :
    FixWallBounceback(_lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | MASK_MASK;
  datamask_modify = X_MASK | V_MASK;

  // copy the blocks (bounds and wall velocities) to the device once

  k_block = DAT::tdual_double_1d("wall/bounceback:block", block.size());
  auto h_block = k_block.view_host();
  for (std::size_t k = 0; k < block.size(); k++) h_block(k) = block[k];
  k_block.modify_host();
  k_block.template sync<DeviceType>();
  d_block = k_block.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   the check for atoms that start inside a wall runs on the host
------------------------------------------------------------------------- */

template <class DeviceType> void FixWallBouncebackKokkos<DeviceType>::setup(int vflag)
{
  atomKK->sync(Host, X_MASK | MASK_MASK);
  FixWallBounceback::setup(vflag);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> void FixWallBouncebackKokkos<DeviceType>::post_integrate()
{
  atomKK->sync(execution_space, datamask_read);
  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  dt = update->dt;

  s_BOUNCES count;
  copymode = 1;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType, TagFixWallBouncebackPostIntegrate>(0, atom->nlocal), *this,
      count);
  copymode = 0;
  atomKK->modified(execution_space, datamask_modify);

  nbounce = count.nbounce;
  nfallback += count.nfallback;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION void FixWallBouncebackKokkos<DeviceType>::operator()(
    TagFixWallBouncebackPostIntegrate, const int &i, s_BOUNCES &count) const
{
  if (!(mask[i] & groupbit)) return;
  double p[3] = {double(x(i, 0)), double(x(i, 1)), double(x(i, 2))};
  double n[3], depth;
  const bool rounded = (mask[i] & slipbit) && rround > 0.0;
  const int b = rounded ? inside_round_kk(p, n, depth) : inside_kk(p);
  if (b < 0) return;

  double xold[3], vold[3], vnew[3];
  int face = -1;
  for (int k = 0; k < 3; k++) {
    vold[k] = double(v(i, k));
    xold[k] = p[k] - dt * vold[k];
  }
  const double s = rounded ? 0.0 : entry_kk(b, xold, p, face);
  if (rounded) {
    // mirror at the rounded surface: reverse the normal velocity only
    double vn = 0.0;
    for (int k = 0; k < 3; k++) vn += (vold[k] - d_block(9 * b + 6 + k)) * n[k];
    for (int k = 0; k < 3; k++) {
      p[k] += 2.0 * depth * n[k];
      vnew[k] = (vn < 0.0) ? vold[k] - 2.0 * vn * n[k] : vold[k];
    }
  } else if ((mask[i] & slipbit) && face >= 0) {
    // mirror at the face: reverse the normal velocity only
    const int k = face % 3;
    for (int m = 0; m < 3; m++) vnew[m] = vold[m];
    p[k] = 2.0 * d_block(9 * b + face) - p[k];
    vnew[k] = 2.0 * d_block(9 * b + 6 + k) - vold[k];
  } else {
    for (int k = 0; k < 3; k++) {
      vnew[k] = 2.0 * d_block(9 * b + 6 + k) - vold[k];
      p[k] = xold[k] + s * (p[k] - xold[k]) + (1.0 - s) * dt * vnew[k];
    }
  }
  if ((rounded ? inside_round_kk(p, n, depth) : inside_kk(p)) >= 0) {
    for (int k = 0; k < 3; k++) {
      p[k] = xold[k];
      vnew[k] = 2.0 * d_block(9 * b + 6 + k) - vold[k];
    }
    count.nfallback += 1.0;
  }
  for (int k = 0; k < 3; k++) {
    x(i, k) = p[k];
    v(i, k) = vnew[k];
  }
  count.nbounce += 1.0;
}

/* ----------------------------------------------------------------------
   index of the first block that contains point p, -1 if none
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION int FixWallBouncebackKokkos<DeviceType>::inside_kk(const double *p) const
{
  for (int b = 0; b < nblock; b++) {
    const int o = 9 * b;
    if (p[0] > d_block(o) && p[0] < d_block(o + 3) && p[1] > d_block(o + 1) &&
        p[1] < d_block(o + 4) && p[2] > d_block(o + 2) && p[2] < d_block(o + 5))
      return b;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   the same with the edges of the blocks rounded with radius rround
   (as in FixWallBounceback::inside_round)
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION int FixWallBouncebackKokkos<DeviceType>::inside_round_kk(const double *p,
                                                                               double *n,
                                                                               double &depth) const
{
  for (int b = 0; b < nblock; b++) {
    const int o = 9 * b;
    double d[3], r2 = 0.0;
    for (int k = 0; k < 3; k++) {
      double q = p[k];
      if (q < d_block(o + k) + rround) q = d_block(o + k) + rround;
      if (q > d_block(o + 3 + k) - rround) q = d_block(o + 3 + k) - rround;
      d[k] = p[k] - q;
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
   fraction of the step xold -> xnew after which the path enters block b,
   and the face it crosses there (as in FixWallBounceback::entry)
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION double FixWallBouncebackKokkos<DeviceType>::entry_kk(int b, const double *xold,
                                                                           const double *xnew,
                                                                           int &face) const
{
  const int o = 9 * b;
  double s = 0.0, t;
  face = -1;
  for (int k = 0; k < 3; k++) {
    if (xold[k] <= d_block(o + k)) {
      t = (d_block(o + k) - xold[k]) / (xnew[k] - xold[k]);
      if (face < 0 || t > s) {
        s = t;
        face = k;
      }
    } else if (xold[k] >= d_block(o + 3 + k)) {
      t = (d_block(o + 3 + k) - xold[k]) / (xnew[k] - xold[k]);
      if (face < 0 || t > s) {
        s = t;
        face = 3 + k;
      }
    }
  }
  if (s < 0.0) s = 0.0;
  return (s < 1.0) ? s : 1.0;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixWallBouncebackKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallBouncebackKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
