/* -*- c++ -*- ----------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.

   If you use this code, please cite Z. Chai, M. de Person, P. A. Buffet,
   M. Franco and G. E. Karniadakis, Biophys. J. (2026),
   https://doi.org/10.1016/j.bpj.2026.10.021.  LAMMPS prints this
   reminder at the end of a run and writes the BibTeX entries to the
   log file (see rbc_cite.h).
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
