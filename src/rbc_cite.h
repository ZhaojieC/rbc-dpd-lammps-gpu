/* -*- c++ -*- ----------------------------------------------------------
   RBC membrane model for LAMMPS
   https://github.com/ZhaojieC/rbc-dpd-lammps
   Distributed under the GNU General Public License, version 2.

   If you use this code, please cite Z. Chai, M. de Person, P. A. Buffet,
   M. Franco and G. E. Karniadakis, Biophys. J. (2026),
   https://doi.org/10.1016/j.bpj.2026.10.021.

   The articles to cite.  Every style of this repository registers them
   with LAMMPS (CiteMe), which prints them at the end of a run and
   writes the BibTeX entries to the log file (log.lammps, or the file
   given with -cite).  Identical entries are listed once, however many
   styles register them.
------------------------------------------------------------------------- */

#ifndef LMP_RBC_CITE_H
#define LMP_RBC_CITE_H

#include "citeme.h"
#include "lammps.h"

namespace LAMMPS_NS {

// The parameters of the model (Table S1) and this implementation.
inline constexpr char CITE_RBC_CHAI2026[] =
    "RBC model of rbc-dpd-lammps: https://doi.org/10.1016/j.bpj.2026.10.021\n\n"
    "@Article{Chai2026,\n"
    " author =  {Z. Chai and M. de Person and P. A. Buffet and M. Franco and G. E. Karniadakis},\n"
    " title =   {Quantifying the biophysical properties of red blood cells in {G}aucher disease},\n"
    " journal = {Biophys.\\ J.},\n"
    " year =    2026,\n"
    " DOI =     {10.1016/j.bpj.2026.10.021}\n"
    "}\n\n";

// The membrane model: spectrin network, area and volume constraints, bending.
inline constexpr char CITE_RBC_FEDOSOV2010[] =
    "RBC membrane model: https://doi.org/10.1016/j.bpj.2010.02.002\n\n"
    "@Article{Fedosov2010,\n"
    " author =  {D. A. Fedosov and B. Caswell and G. E. Karniadakis},\n"
    " title =   {A multiscale red blood cell model with accurate mechanics, rheology, and dynamics},\n"
    " journal = {Biophys.\\ J.},\n"
    " year =    2010,\n"
    " volume =  98,\n"
    " number =  10,\n"
    " pages =   {2215--2225},\n"
    " DOI =     {10.1016/j.bpj.2010.02.002}\n"
    "}\n\n";

// Register the articles to cite; the membrane styles also cite the model.
inline void rbc_cite(LAMMPS *lmp, bool membrane)
{
  if (!lmp->citeme) return;
  lmp->citeme->add(CITE_RBC_CHAI2026);
  if (membrane) lmp->citeme->add(CITE_RBC_FEDOSOV2010);
}

}    // namespace LAMMPS_NS

#endif
