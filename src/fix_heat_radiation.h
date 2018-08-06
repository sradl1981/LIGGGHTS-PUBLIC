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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2016-     Thomas Forgber, TU Graz
------------------------------------------------------------------------- */
#ifdef FIX_CLASS

FixStyle(heat/radiation,FixHeatRad)

#else

#ifndef LMP_FIX_HEATRAD_H
#define LMP_FIX_HEATRAD_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixHeatRad : public Fix {

  public:

    FixHeatRad(class LAMMPS *, int, char **);

    ~FixHeatRad();

    virtual void post_create();
    virtual void init();
    virtual int setmask();
    void updatePtrs();  
    void initial_integrate(int vflag);
    virtual void final_integrate();
    double FindDistanceToSegment(double x1, double y1, double z1, double x2, double y2, double z2, double pointX, double pointY, double pointZ);
    double t;

  /*  virtual void pre_delete(bool unfixflag){ UNUSED(unfixflag); };
  
    virtual double compute_scalar();
    
    // per default these three methods throw errors.
    virtual void cpl_evaluate(class ComputePairGranLocal *);
    virtual void register_compute_pair_local(class ComputePairGranLocal *);
    virtual void unregister_compute_pair_local(class ComputePairGranLocal *);
    */
  private:

    //vector< vector<double> > view_factors_;

    bool init_points_on_sphere_;

    double* x_points_origin_;
    double* y_points_origin_;
    double* z_points_origin_;
    int number_of_distributed_points_;

    int  solve_every_ ;
    bool solve_reflection_; 
    int  solve_counter_;    
    /*double x_points_particle_[];
    double y_points_particle_[];
    double z_points_particle_[]; */   
    
  protected:

    class ComputePairGranLocal *cpl;
    class FixPropertyAtom* fix_heatFlux;
    class FixPropertyAtom* fix_heatSource;
    class FixPropertyAtom* fix_temp;
    class FixScalarTransportEquation *fix_ste;
    class FixPropertyAtom* fix_raddirectionalHeatFlux;
	class FixPropertyAtom* fix_accumulativeviewfactor;
    class PairGran *pair_gran;

    double *heatFlux;   
    int history_flag;
    int number_of_points_;
    double emissivity_;
	double radius_drum_;
	double temp_drum_;
    double *heatSource; 
    double *Temp;       
    double T0;          
    double **raddirectionalHeatFlux;
	double *accumulativeviewfactor;
     
  };

}

#endif
#endif
