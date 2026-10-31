/* -*- c++ -*- ----------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(area/volume,AngleAreaVolume);
// clang-format on
#else

#ifndef LMP_ANGLE_AREA_VOLUME_H
#define LMP_ANGLE_AREA_VOLUME_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleAreaVolume : public Angle {
 public:
  AngleAreaVolume(class LAMMPS *);
  ~AngleAreaVolume() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_angle(int) override { return -1.0; }    // not an angle potential
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, int, int, int) override;
  void *extract(const char *, int &) override;
  double memory_usage() override;

  // area and volume of every cell from the current positions, indexed by
  // molecule ID; returns the largest molecule ID
  int cell_totals(double *&, double *&, int *&);

 protected:
  double *ka, *area0, *kv, *volume0, *kd;
  class FixStoreAtom *ref;

  // one entry per molecule ID 0..maxmol (a cell is one molecule)
  int maxmol;
  int *celltype;          // angle type used by the molecule, 0 if it has no angles
  double *nmolatoms;      // number of atoms in the molecule
  double *com;            // reference point of the molecule (3 per molecule)
  double *area, *volume;  // current total area and enclosed volume
  double *buf1, *buf2;    // reduction buffers

  void allocate();
  void setup_cells();
  void reference_points();
  void totals();
  double triangle_area0(int, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
