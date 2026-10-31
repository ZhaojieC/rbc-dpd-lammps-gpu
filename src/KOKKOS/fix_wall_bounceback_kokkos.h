/* -*- c++ -*- ----------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps-gpu
   Distributed under the GNU General Public License, version 2.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/bounceback/kk,FixWallBouncebackKokkos<LMPDeviceType>);
FixStyle(wall/bounceback/kk/device,FixWallBouncebackKokkos<LMPDeviceType>);
FixStyle(wall/bounceback/kk/host,FixWallBouncebackKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_WALL_BOUNCEBACK_KOKKOS_H
#define LMP_FIX_WALL_BOUNCEBACK_KOKKOS_H

#include "fix_wall_bounceback.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixWallBouncebackPostIntegrate {};

// reflections in one step, summed over atoms
struct s_BOUNCES {
  double nbounce, nfallback;
  KOKKOS_INLINE_FUNCTION
  s_BOUNCES() : nbounce(0.0), nfallback(0.0) {}
  KOKKOS_INLINE_FUNCTION
  s_BOUNCES &operator+=(const s_BOUNCES &rhs)
  {
    nbounce += rhs.nbounce;
    nfallback += rhs.nfallback;
    return *this;
  }
};

template <class DeviceType> class FixWallBouncebackKokkos : public FixWallBounceback {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixWallBouncebackKokkos(class LAMMPS *, int, char **);
  void setup(int) override;
  void post_integrate() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallBouncebackPostIntegrate, const int &, s_BOUNCES &) const;

 protected:
  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_int_1d_randomread mask;
  DAT::tdual_double_1d k_block;
  typename AT::t_double_1d d_block;
  double dt;

  KOKKOS_INLINE_FUNCTION
  int inside_kk(const double *) const;
  KOKKOS_INLINE_FUNCTION
  int inside_round_kk(const double *, double *, double &) const;
  KOKKOS_INLINE_FUNCTION
  double entry_kk(int, const double *, const double *, int &) const;
};

}    // namespace LAMMPS_NS

#endif
#endif
