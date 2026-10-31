/* -*- c++ -*- ----------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(wlc/pow,BondWLCPOW);
// clang-format on
#else

#ifndef LMP_BOND_WLC_POW_H
#define LMP_BOND_WLC_POW_H

#include "bond.h"

namespace LAMMPS_NS {

class BondWLCPOW : public Bond {
 public:
  BondWLCPOW(class LAMMPS *);
  ~BondWLCPOW() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_distance(int) override { return -1.0; }    // depends on the bond
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *ratio;    // lmax / l0
  double *mu0;      // target shear modulus of the network
  double *expo;     // exponent m of the repulsive POW term
  class FixStoreAtom *ref;

  void allocate();
  double reference_length(int, int);
  void spring(int, double, double &, double &, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif
