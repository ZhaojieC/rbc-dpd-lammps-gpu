/* -*- c++ -*- ----------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(rbc/shape,ComputeRBCShape);
// clang-format on
#else

#ifndef LMP_COMPUTE_RBC_SHAPE_H
#define LMP_COMPUTE_RBC_SHAPE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRBCShape : public Compute {
 public:
  ComputeRBCShape(class LAMMPS *, int, char **);
  ~ComputeRBCShape() override;
  void init() override;
  void compute_array() override;

 private:
  int nrows;
  double **data;
  class AngleAreaVolume *av;
};

}    // namespace LAMMPS_NS

#endif
#endif
