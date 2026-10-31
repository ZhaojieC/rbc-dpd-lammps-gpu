/* -*- c++ -*- ----------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.
------------------------------------------------------------------------- */

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(bend,DihedralBend);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_BEND_H
#define LMP_DIHEDRAL_BEND_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralBend : public Dihedral {
 public:
  DihedralBend(class LAMMPS *);
  ~DihedralBend() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  double *k, *theta0, *cos0, *sin0;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
