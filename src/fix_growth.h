/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Copyright 2017      Stefan Radl, TU Graz, 2017
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(growth,FixGrowth)

#else

#ifndef LMP_FIX_GROWTH_H
#define LMP_FIX_GROWTH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGrowth : public Fix {
 public:

  FixGrowth(class LAMMPS *, int, char **);
  ~FixGrowth();
  int setmask();
  void post_create(); 
  void pre_delete(bool unfixflag); 

  void setup_pre_force(int);
  void pre_force(int);

 private:

  class FixAdapt *fixAdapt_;
  int  ifixAdapt_;
  char fixid[100];
  char *nEveryString_;
  char *variableToControlGrowth_;

  double updateInterval_;

  bool   doDynamicGrowth_;

  void   change_variableForGrowth();
  void   initDynamicGrowth();

  //Fixes for dynamic growth modeling
  char *fix_name_supersaturation_;
  char *fix_name_temperature_;
  char *fix_name_sherwood_;
  char *fix_name_saturationdensity_;
  class FixPropertyAtom *fix_supersaturation_;
  class FixPropertyAtom *fix_temperature_;
  class FixPropertyAtom *fix_sherwood_;
  class FixPropertyAtom *fix_saturationdensity_;
  class FixPropertyGlobal* fix_diffusionConstant_;
  class FixPropertyGlobal* fix_surfaceTensionConstant_;

};

}

#endif
#endif

