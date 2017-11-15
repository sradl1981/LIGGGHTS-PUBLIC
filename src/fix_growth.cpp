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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_growth.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h" 
#include "modify.h" 
#include "update.h" 
#include "fix_adapt.h" 

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixGrowth::FixGrowth(LAMMPS *lmp, int narg, char **arg) : 
Fix(lmp, narg, arg)
{
    if (narg < 5) error->all(FLERR,"Illegal fix adapt command");
    nevery = force->inumeric(FLERR,arg[3]);
    if (nevery < 0) error->all(FLERR,"Illegal fix adapt command");
 
    //save the nevery as a string for future use
    int nA = strlen(&arg[3][0]) + 1;
    nEveryString_ = new char[nA];
    strcpy(nEveryString_,&arg[3][0]);

    int iarg = 4;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix adapt command");
      if (strcmp(arg[iarg+1],"diameter") == 0) {
        //XXX: do nothing now, use in future to make settings if necessary
      } else error->all(FLERR,"Illegal fix adapt command");
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][0]) + 1;
        variableToControlGrowth_ = new char[n];
        strcpy(variableToControlGrowth_,&arg[iarg+2][0]);
      } else error->all(FLERR,"Illegal fix adapt command");
      iarg += 3;
    } else break;
  }

  //generate the growth variable
  char *varargs[3];
  varargs[0]=&variableToControlGrowth_[2];
  varargs[1]=(char *)"atom";
  varargs[2]=(char *)"1e-4";

  input->variable->set(3,varargs); //add variable to global set of vars
  printf("FixGrowth generated the variable with name: %s. \n", varargs[0]);
  fixAdapt_ = NULL; 
}

/* ---------------------------------------------------------------------- */

FixGrowth::~FixGrowth()
{

}

/* ---------------------------------------------------------------------- */

void FixGrowth::post_create()
{

  if (fixAdapt_ == NULL)
  {
    const char *fixarg[7];
    
    sprintf(fixid,"growth_%s",id);
    fixarg[0]=fixid;
    fixarg[1]="all";
    fixarg[2]="adapt";
    fixarg[3]=nEveryString_;
    fixarg[4]="atom"; 
    fixarg[5]="diameter";    
    fixarg[6]=variableToControlGrowth_;  
    modify->add_fix(7,const_cast<char**>(fixarg));

      int ifix = modify->find_fix(fixid);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for 'adapt' fix does not exist in the growth fix.");

  }
}

/* ---------------------------------------------------------------------- */

void FixGrowth::pre_delete(bool unfixflag)
{
    if (unfixflag && fixAdapt_) modify->delete_fix(fixid); 
}

/* ---------------------------------------------------------------------- */

int FixGrowth::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_RUN;
  return mask;
}



