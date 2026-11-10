/* ----------------------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.

   If you use this code, please cite Z. Chai, M. de Person, P. A. Buffet,
   M. Franco and G. E. Karniadakis, Biophys. J. (2026),
   https://doi.org/10.1016/j.bpj.2026.10.021.  LAMMPS prints this
   reminder at the end of a run and writes the BibTeX entries to the
   log file (see rbc_cite.h).
------------------------------------------------------------------------- */

#include "rbc_reference.h"

#include "accelerator_kokkos.h"
#include "atom.h"
#include "atom_masks.h"
#include "error.h"
#include "fix_store_atom.h"
#include "lammps.h"
#include "modify.h"

using namespace LAMMPS_NS;

static const char *const REF_ID = "RBC_X0";

FixStoreAtom *LAMMPS_NS::rbc_reference(LAMMPS *lmp)
{
  auto *ref = dynamic_cast<FixStoreAtom *>(lmp->modify->get_fix_by_id(REF_ID));
  if (!ref) {
    // after read_restart, add_fix() fills the new fix from the restart file
    ref = dynamic_cast<FixStoreAtom *>(
        lmp->modify->add_fix(std::string(REF_ID) + " all STORE/ATOM 4 0 1 1"));
    if (!ref) lmp->error->all(FLERR, "Could not create the RBC reference shape storage");
    ref->create_attribute = 1;    // atoms created later start with flag 0
  }

  // with KOKKOS the host copy of the coordinates may be out of date
  if (lmp->kokkos) dynamic_cast<AtomKokkos *>(lmp->atom)->sync(Host, X_MASK);

  double **x0 = ref->astore;
  double **x = lmp->atom->x;
  const int nlocal = lmp->atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (x0[i][3] == 0.0) {
      x0[i][0] = x[i][0];
      x0[i][1] = x[i][1];
      x0[i][2] = x[i][2];
      x0[i][3] = 1.0;
    }
  }
  return ref;
}
