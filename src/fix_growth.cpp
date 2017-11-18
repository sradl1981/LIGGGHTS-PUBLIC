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
    doDynamicGrowth_ = false;

    if (narg < 5) error->all(FLERR,"Illegal fix adapt command");
    nevery = force->inumeric(FLERR,arg[3]);
    if (nevery < 0) error->all(FLERR,"Illegal fix adapt command");
 
    //save the nevery as a string for future use, as well as the updateInterval_
    int nA = strlen(&arg[3][0]) + 1;
    nEveryString_ = new char[nA];
    strcpy(nEveryString_,&arg[3][0]);
    updateInterval_ = update->dt * atof(nEveryString_);
    printf("FixGrowth will use the following updateInterval: %g. \n", updateInterval_);

    int iarg = 4;
    char *growToValue;
    growToValue = new char[1]; growToValue = (char *)"0";

    while (iarg < narg) {
        if (strcmp(arg[iarg],"atom") == 0) {
              if (iarg+3 > narg) error->all(FLERR,"Illegal fix adapt command");
              if (strcmp(arg[iarg+1],"diameter") == 0) {
                //do nothing now, may use in future to make other settings if necessary
              } else error->all(FLERR,"Illegal fix adapt command");
              if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
                int n = strlen(&arg[iarg+2][0]) + 1;
                variableToControlGrowth_ = new char[n];
                strcpy(variableToControlGrowth_,&arg[iarg+2][0]);
              } else error->all(FLERR,"Illegal fix adapt command. Cannot find 'v_' keyword.");
              if (strstr(arg[iarg+3],"growToValue") == arg[iarg+3]) {
                if(narg<9) error->all(FLERR,"Input after keyword 'growToValue' expected .");
                int n = strlen(&arg[iarg+4][0]) + 1;
                growToValue = new char[n];
                strcpy(growToValue,&arg[iarg+4][0]);
              } else error->all(FLERR,"Illegal fix adapt command. Cannot find 'growToValue' keyword.");
        iarg += 3;
    } else break;
  }

  //generate a compute that calculates the particle radius
  char *computearg[4];
  computearg[0]=(char *)"radius";
  computearg[1]=(char *)"all";
  computearg[2]=(char *)"property/atom";
  computearg[3]=(char *)"radius";
  modify->add_compute(4,computearg,lmp->suffix);

  //generate the growth variable, and activate dynamic growth if necessary
  fix_name_supersaturation_ =  new char[20];
  fix_name_temperature_ =  new char[20];
  fix_name_sherwood_ =  new char[20];
  fix_name_saturationdensity_ =  new char[20];

  char *varargs[3];
  varargs[0]=&variableToControlGrowth_[2];
  varargs[1]=(char *)"atom";
  if(strcmp(growToValue,"DYNAMIC")==0) 
  { 
    strcpy(fix_name_supersaturation_, "supersat"); //TODO:hardcoded, might want to read in
    strcpy(fix_name_temperature_, "temp"); //TODO:hardcoded, might want to read in
    strcpy(fix_name_sherwood_, "sherwood"); //TODO:hardcoded, might want to read in
    strcpy(fix_name_saturationdensity_, "saturationdensity"); //TODO:hardcoded, might want to read in
    initDynamicGrowth();
    varargs[2]=(char *)"2.0*c_radius"; //do not grow, just keep the same for now
  }
  else
  {
      strcpy(fix_name_supersaturation_, "N/A");
      strcpy(fix_name_temperature_, "N/A");
      strcpy(fix_name_sherwood_, "N/A");
      strcpy(fix_name_saturationdensity_, "N/A");
      varargs[2]=growToValue;
  }
  input->variable->set(3,varargs); //add variable to global set of vars
  printf("FixGrowth generated the variable with name: %s. \n", varargs[0]);
  fixAdapt_ = NULL; 
  ifixAdapt_= -1;

  fix_supersaturation_ = NULL;
  fix_temperature_     = NULL;
  fix_diffusionConstant_ = NULL;
  fix_surfaceTensionConstant_ = NULL;

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
    
    sprintf(fixid,"fixGrowth_%s",id);
    fixarg[0]=fixid;
    fixarg[1]="all";
    fixarg[2]="adapt";
    fixarg[3]=nEveryString_;
    fixarg[4]="atom"; 
    fixarg[5]="diameter";  
    fixarg[6]=variableToControlGrowth_;  
    modify->add_fix(7,const_cast<char**>(fixarg));

      ifixAdapt_ = modify->find_fix(fixid);
      if (ifixAdapt_ < 0)
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

/* ---------------------------------------------------------------------- */
void FixGrowth::setup_pre_force(int vflag)
{
  change_variableForGrowth();
}

/* ---------------------------------------------------------------------- */
void FixGrowth::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;
  change_variableForGrowth();
}

/* ---------------------------------------------------------------------- */
void FixGrowth::change_variableForGrowth()
{
    if(!doDynamicGrowth_) return;

    printf("FixGrowth::change_variableForGrowth(). \n");
    //This implements the growth physics!
    //must loop over the particles, and set the variable for growth individually
}

/* ---------------------------------------------------------------------- */
void FixGrowth::initDynamicGrowth()
{
    doDynamicGrowth_ = true;

    //ensure all fixes that specify dynamic growth are available
    fix_supersaturation_ =
        static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_name_supersaturation_,"property/atom","scalar",0,0,style));

    fix_temperature_ =
        static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_name_temperature_,"property/atom","scalar",0,0,style));

    fix_sherwood_ =
        static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_name_sherwood_,"property/atom","scalar",0,0,style));

    fix_saturationdensity_ =
        static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_name_saturationdensity_,"property/atom","scalar",0,0,style));

    fix_diffusionConstant_=static_cast<FixPropertyGlobal*>(modify->find_fix_property("diffusionCoeff","property/global","scalar",0,0,style));

    fix_surfaceTensionConstant_=static_cast<FixPropertyGlobal*>(modify->find_fix_property("surfaceTension","property/global","scalar",0,0,style));

}


