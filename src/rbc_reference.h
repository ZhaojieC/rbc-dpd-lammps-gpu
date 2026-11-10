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

#ifndef LMP_RBC_REFERENCE_H
#define LMP_RBC_REFERENCE_H

namespace LAMMPS_NS {

class LAMMPS;
class FixStoreAtom;

// Reference (stress-free) membrane shape shared by bond wlc/pow and
// angle area/volume.  The internal fix RBC_X0 stores, for every atom,
// its reference position (columns 0-2) and a flag (column 3) that is
// set once the position has been recorded.  Atoms that have no
// reference yet take their current position, so the reference shape is
// the shape at the first run after the atoms were read or created.
// The fix writes its data to restart files and read_restart restores it.

FixStoreAtom *rbc_reference(LAMMPS *lmp);

}    // namespace LAMMPS_NS

#endif
