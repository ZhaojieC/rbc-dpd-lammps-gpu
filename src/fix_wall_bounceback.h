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

#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/bounceback,FixWallBounceback);
// clang-format on
#else

#ifndef LMP_FIX_WALL_BOUNCEBACK_H
#define LMP_FIX_WALL_BOUNCEBACK_H

#include "fix.h"

#include <vector>

namespace LAMMPS_NS {

class FixWallBounceback : public Fix {
 public:
  FixWallBounceback(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_integrate() override;
  double compute_vector(int) override;

 protected:
  int nblock;
  std::vector<double> block;    // per block: lo[3], hi[3], wall velocity[3]
  int slipbit;                  // group bit of the particles reflected like a mirror
  double rround;                // radius of the block edges seen by the slip group
  double nbounce;               // reflections on this processor in the last step
  double nfallback;             // reflections that had to restore the old position

  int inside(const double *) const;
  int inside_round(const double *, double *, double &) const;
  int inwall(const double *, int) const;
  double entry(int, const double *, const double *, int &) const;
};

}    // namespace LAMMPS_NS

#endif
#endif
