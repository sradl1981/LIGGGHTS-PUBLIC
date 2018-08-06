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

#include "fix_heat_radiation.h"

#include <math.h>
#include "atom.h"
#include "fix_property_atom.h"
#include "fix_scalar_transport_equation.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "modify.h"
#include "pair_gran.h"
#include <stdlib.h>
#include "neigh_list.h"
#include "fix_heat_gran_conduction.h"
#include "compute_pair_gran_local.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "math_extra_liggghts.h"
#include "math_const.h"
#include "properties.h"
#include "mpi.h"
#include "domain.h"

#define SMALL_VIEW_FACTOR 1.e-16
#define SMALL_NUMBER 1.e-16
#define LIMIT_DRUM_RADIUS 10   		//hardcode limit for radius of drum


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixHeatRad::FixHeatRad(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  init_points_on_sphere_(false),
  number_of_distributed_points_(0),
  solve_every_(0),
  solve_reflection_(false),
  solve_counter_(0),
  radius_drum_(0.0),
  temp_drum_(0.0)
{
 //pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix heat/gran needs per particle radius and mass");

  if (narg < 5)
    error->fix_error(FLERR,this,"not enough arguments");

  int iarg = 3;

  if(strcmp(arg[iarg++],"number_of_points"))
    error->fix_error(FLERR,this,"expecting keyword 'number_of_points'");
  number_of_points_ = atof(arg[iarg++]);
 
  if(strcmp(arg[iarg++],"emissivity"))
    error->fix_error(FLERR,this,"expecting keyword 'emissivity'");
  emissivity_ = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"rotating_drum_radius"))
    error->fix_error(FLERR,this,"expecting keyword 'rotating_drum_radius'");
    radius_drum_ = atof(arg[iarg++]); 

  if(strcmp(arg[iarg++],"rotating_drum_temperature"))
    error->fix_error(FLERR,this,"expecting keyword 'rotating_drum_temperature'");
    temp_drum_ = atof(arg[iarg++]); 

  if(strcmp(arg[iarg++],"solve_every"))
    error->fix_error(FLERR,this,"expecting keyword 'solve_every'");
  solve_every_ = atoi(arg[iarg++]);

  //if(strcmp(arg[iarg++],"solvereflection"))
  //  solve_reflection_ = true;
    //error->fix_error(FLERR,this,"expecting keyword 'emissivity'");
   
  fix_temp = fix_heatFlux = fix_heatSource = NULL;
  fix_ste = NULL;
  fix_raddirectionalHeatFlux = NULL;
  fix_accumulativeviewfactor = NULL;

}

/* ---------------------------------------------------------------------- */

FixHeatRad::~FixHeatRad()
{
    //view_factors_.clear();
   //Free points and particles 
    delete [] x_points_origin_;
    delete [] y_points_origin_;
    delete [] z_points_origin_;
}

/* ---------------------------------------------------------------------- */

void FixHeatRad::post_create()
{
  // register directional flux
  fix_raddirectionalHeatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("raddirectionalHeatFlux","property/atom","vector",3,0,this->style,false));
  if(!fix_raddirectionalHeatFlux)
  {
    const char* fixarg[11];
    fixarg[0]="raddirectionalHeatFlux";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="raddirectionalHeatFlux";
    fixarg[4]="vector";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fix_raddirectionalHeatFlux = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
  }

  // register accumulative view factor
  fix_accumulativeviewfactor = static_cast<FixPropertyAtom*>(modify->find_fix_property("accumulativeviewfactor","property/atom","scalar",0,0,this->style,false));
  if(!fix_accumulativeviewfactor)
  {
    const char* fixarg[9];
    fixarg[0]="accumulativeviewfactor";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="accumulativeviewfactor";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fix_accumulativeviewfactor = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatRad::updatePtrs(){
	
  Temp = fix_temp->vector_atom;
  vector_atom = Temp; 

  heatFlux = fix_heatFlux->vector_atom;
  heatSource = fix_heatSource->vector_atom;
  raddirectionalHeatFlux = fix_raddirectionalHeatFlux->array_atom;
  accumulativeviewfactor = fix_accumulativeviewfactor->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixHeatRad::init()
{

  if (!atom->radius_flag || !atom->rmass_flag)
    error->fix_error(FLERR,this,"must use a granular atom style ");

    // check if a fix of this style already exists
  if(modify->n_fixes_style(style) > 1)
    error->fix_error(FLERR,this,"cannot have more than one fix of this style");

  if(!force->pair_match("gran", 0))
    error->fix_error(FLERR,this,"needs a granular pair style to be used");

   pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
   history_flag = pair_gran->is_history();

  fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste) error->fix_error(FLERR,this,"Needs a fix transportequation/scalar to perform the integration. Please ensure such a fix is available by loading a fix of type 'heat/gran/conduction'");

  fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
  fix_heatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));
  fix_heatSource = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatSource","property/atom","scalar",0,0,style));
  fix_raddirectionalHeatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("raddirectionalHeatFlux","property/atom","vector",0,0,style));
  fix_accumulativeviewfactor = static_cast<FixPropertyAtom*>(modify->find_fix_property("accumulativeviewfactor","property/atom","scalar",0,0,style));

}

/* ---------------------------------------------------------------------- */

int FixHeatRad::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatRad::initial_integrate(int vflag)
{
    if(solve_counter_ == 0)
    {
      solve_counter_ = solve_every_;
      double dirFlux[3];
      double flux;
      int i,j,ii,jj,inum,jnum;
      double x1tmp,y1tmp,z1tmp;
      double x2tmp,y2tmp,z2tmp;
      double rad1i,rad2j=0.0;
      int *ilist,*jlist,*numneigh,**firstneigh;
      double delx, dely, delz;
    
      double x_points_particle_[number_of_points_+2];
      double y_points_particle_[number_of_points_+2];
      double z_points_particle_[number_of_points_+2];
      
        inum = pair_gran->list->inum;
        ilist = pair_gran->list->ilist;
        numneigh = pair_gran->list->numneigh;
        firstneigh = pair_gran->list->firstneigh;

       double *radius = atom->radius;
       double **x = atom->x;
       double a_, d_ , d_theta_, d_phi_, theta_,  phi_ ;
       double radius_mono_ = radius[ilist[0]];

       int M_phi_ , M_theta_ ;
       double M_;
       double S_particle_ ;

       double correction_a_ =  0.03407;
       double correction_b_ = -1.965;
       double distance_c_c_; 
       double correction_angle_;
       double solid_angle_;
       double maximum_distance_from_c_c_vec_;
       updatePtrs();

       //Point distribution - just do once and update positions with particle positions
        if(!init_points_on_sphere_)
        {
           x_points_origin_ = new double[number_of_points_+5];
           y_points_origin_ = new double[number_of_points_+5];
           z_points_origin_ = new double[number_of_points_+5];

           a_ = MathConst::MY_4PI * radius_mono_ * radius_mono_ /  number_of_points_ ;
           d_ = sqrt(a_);
           M_theta_ = round(MathConst::MY_PI * radius_mono_ / d_) ;
           d_theta_ = radius_mono_ *  MathConst::MY_PI / M_theta_;
           d_phi_ = a_ / d_theta_ ; 
            
           for (int m_ = 0; m_ < (M_theta_) ; m_++) 
           {
              M_ = m_+0.5;
              theta_  = MathConst::MY_PI*(M_)/M_theta_ ;
              M_phi_ =  round(radius_mono_*MathConst::MY_2PI*sin(theta_)/d_phi_);
                
              for (int n_ = 0; n_ < M_phi_ ; n_++) 
              {
                  phi_ =  MathConst::MY_2PI * n_ / M_phi_;

                  if(number_of_distributed_points_ >= (number_of_points_+4))
                    error->fix_error(FLERR,this,"You try to distribute more points than you assigned");

                  x_points_origin_[number_of_distributed_points_] = radius_mono_ * (sin(theta_)*cos(phi_));
                  y_points_origin_[number_of_distributed_points_] = radius_mono_ * (sin(theta_)*sin(phi_));
                  z_points_origin_[number_of_distributed_points_] = radius_mono_ * (cos(theta_));
                  number_of_distributed_points_++;
                  
              }
           }
          init_points_on_sphere_ = true;
       }   

       for (ii = 0; ii < inum; ii++) 
       {
            raddirectionalHeatFlux[ii][0] = 0.0;
            raddirectionalHeatFlux[ii][1] = 0.0;
            raddirectionalHeatFlux[ii][2] = 0.0;
			accumulativeviewfactor[ii] = 0.0;
       }
       fix_raddirectionalHeatFlux->do_forward_comm();
	   fix_accumulativeviewfactor->do_forward_comm();

       updatePtrs();

       // Loop over all particles 
        for (ii = 0; ii < inum; ii++) 
        {
            i = ilist[ii];
            x1tmp = x[i][0];
            y1tmp = x[i][1];
            z1tmp = x[i][2];
            
            if(Temp[i] < SMALL_NUMBER)
                    Temp[i] = 0.0;

            int view_factor_count_ [number_of_distributed_points_];
            for (int k_ = 0; k_ < number_of_distributed_points_; k_++) 
            {
                //TODO: ensure that  k_ is inside pre-allocated arrays!
                //k_ has to be inside because check for number_of_distributed_points_ is added
                x_points_particle_[k_] = x_points_origin_[k_] + x1tmp;
                y_points_particle_[k_] = y_points_origin_[k_] + y1tmp;
                z_points_particle_[k_] = z_points_origin_[k_] + z1tmp;
                view_factor_count_[k_] = 0;
            }
              
            rad1i = radius[i];
            jlist = firstneigh[i];
            jnum = numneigh[i];
            double view_factor_calc_ [jnum];
            double distances_ [jnum];

            //Loop over neighbor list
            for (jj = 0; jj < jnum; jj++) 
            {
                  j = jlist[jj];
                  j &= NEIGHMASK;
                  view_factor_calc_ [jj] = 0.0;
                  x2tmp = x[j][0];
                  y2tmp = x[j][1];
                  z2tmp = x[j][2];

                  
                  if(domain->is_periodic_ghost(j))
                      //printf("here because particle (%g %g %g) and particle (%g %g %g) are neighbors and inside the domain) \n",x1tmp,y1tmp,z1tmp,x2tmp,y2tmp,z2tmp);
                  {
                      rad2j = radius[j];

                      // calculate distance and correction factor
                      distance_c_c_ = sqrt(  (x2tmp-x1tmp)*(x2tmp-x1tmp) 
                                           + (y2tmp-y1tmp)*(y2tmp-y1tmp)  
                                           + (z2tmp-z1tmp)*(z2tmp-z1tmp)
                                          ) ;

                      //Set min distance in case of overlap since view factor is not defined for penetrating particles
                      if(distance_c_c_ < (rad2j+rad1i))
                      {
                      	distance_c_c_ = rad2j+rad1i;
                      }
                        
                      distances_ [jj] = distance_c_c_;

                      // S_particle_ is dimensionless distance between surfaces - needed for correction
                      S_particle_ = (distance_c_c_ - (rad2j+rad1i)) / rad1i;
        
                      // Calculate correction factor - may be speeded up
                      correction_angle_ = correction_a_*exp(correction_b_*S_particle_); 
					  //correction_angle_ = 0.0;
                      // Calculate solid angle including correction factor 
                      solid_angle_ = asin(rad2j/distance_c_c_)+correction_angle_;
                      maximum_distance_from_c_c_vec_ = sin(solid_angle_)*rad1i;

                      for (int k_ = 0; k_ < number_of_distributed_points_; k_++) 
                      { 
                        // Faster with a large number of particles    
                        double distance_ = FindDistanceToSegment(x1tmp, y1tmp, z1tmp, x2tmp,y2tmp,z2tmp ,x_points_particle_[k_], y_points_particle_[k_], z_points_particle_[k_]);
     
                       // time_ = MPI_Wtime();
                        if(distance_ <= maximum_distance_from_c_c_vec_)
                        {
                           if(view_factor_count_ [k_] == 0 )
                           {
                                view_factor_count_ [k_] = jj+1;                        // Set view factor
                                view_factor_calc_ [jj] = view_factor_calc_ [jj] + 1.0;  // Register point 
                           
                           }
                           else if(distance_c_c_ <= distances_[view_factor_count_ [k_]-1])
                           {
                                view_factor_calc_ [(view_factor_count_ [k_]-1)] = view_factor_calc_ [(view_factor_count_ [k_])-1] - 1.0;  // Unregister point for old particle 
                                view_factor_count_ [k_] = jj+1; 
                                view_factor_calc_ [jj] = view_factor_calc_ [jj] + 1.0;  // Register point for new particle 
                           }
                        }
                    
                      } // end point loop  
                    }
            } // end neightbor loop
        
            for (jj = 0; jj < jnum; jj++) 
            {

                j = jlist[jj];		//Real index of neightbor

                if(Temp[j] < SMALL_NUMBER)
                    Temp[j] = 0.0;

                x2tmp = x[j][0];
                y2tmp = x[j][1];
                z2tmp = x[j][2];

                if(domain->is_periodic_ghost(j))
                {
					//printf("Here because periodic \n");
                    delx = x1tmp - x2tmp;
                    dely = y1tmp - y2tmp;
                    delz = z1tmp - z2tmp;

                    view_factor_calc_ [jj] = view_factor_calc_ [jj]/number_of_distributed_points_;
					accumulativeviewfactor [i] += view_factor_calc_ [jj];
					accumulativeviewfactor [j] += view_factor_calc_ [jj];
					
					
					//if(jj == (jnum-1))
					//{
					//	printf("particle i  = %i , neigh = %i, has accum_view_factor_ [jj] = %g \n",i,jj,accum_view_factor_ [i]);
					//}
                    if(i==0)
                    {
                        //printf("View factor of %i to %i = %g \n", i,j,view_factor_calc_ [jj]);
                    }

                    flux = solve_every_ * (MathConst::MY_BOLZ_CONST*(  Temp[i]*Temp[i]*Temp[i]*Temp[i] 
                                                      - Temp[j]*Temp[j]*Temp[j]*Temp[j])
                                                     )
                            / (
                                 ((1.0-emissivity_)/(MathConst::MY_PI*rad1i*rad1i))
                               + ((1.0-emissivity_)/(MathConst::MY_PI*rad2j*rad2j)) 
                               + (1.0/(4.0*MathConst::MY_PI*rad1i*rad1i*(view_factor_calc_ [jj]+SMALL_VIEW_FACTOR)))
                              );

					dirFlux[0] = flux*delx;
                    dirFlux[1] = flux*dely;
                    dirFlux[2] = flux*delz;

                    raddirectionalHeatFlux[i][0] -= 0.50 * dirFlux[0];		//TODO: CHECK!!!
                    raddirectionalHeatFlux[i][1] -= 0.50 * dirFlux[1];
                    raddirectionalHeatFlux[i][2] -= 0.50 * dirFlux[2];

                    raddirectionalHeatFlux[j][0] += 0.50 * dirFlux[0];
                    raddirectionalHeatFlux[j][1] += 0.50 * dirFlux[1];
                    raddirectionalHeatFlux[j][2] += 0.50 * dirFlux[2];  
	
					

                    //printf("Flux = %g \n",flux);
                    heatFlux[i] -= flux;
                    heatFlux[j] += flux;
                }

            }//end neightbor loop
            
			// Start wall radiation exchange
			if(radius_drum_<LIMIT_DRUM_RADIUS)
					{
						double distance_from_center_ = sqrt(x1tmp*x1tmp + y1tmp*y1tmp);
                        //printf("distance from center = %g \n",distance_from_center_);

						if(distance_from_center_>radius_drum_)
							error->all(FLERR,"Particles are placed outside of the drum");
						double distance_from_wall_ = (radius_drum_-rad1i)-distance_from_center_;
						if(distance_from_wall_<0.0)
							distance_from_wall_=0.0;

						//double pow_dimless_distance_ = (distance_from_wall_/rad1i)*(distance_from_wall_/rad1i);
						
						double view_factor_wall_ = exp(-2.6-0.8*((distance_from_wall_/rad1i)*(distance_from_wall_/rad1i))) 
                                                   + 0.12*exp(-12.0*(distance_from_wall_/rad1i)*(distance_from_wall_/rad1i));

						double view_factor_rest_ = 1.0 - (accumulativeviewfactor[i]+view_factor_wall_);

						
						if(view_factor_rest_<0.0)
						{
							view_factor_wall_=0.0;
						}
	
						printf("Part %i has dimless distance from wall %g, radius %g \n", i,distance_from_wall_/rad1i,rad1i);
						accumulativeviewfactor[i] += view_factor_wall_;

						double flux_wall_ = solve_every_ * (MathConst::MY_BOLZ_CONST*(  Temp[i]*Temp[i]*Temp[i]*Temp[i] 
                                                      - temp_drum_*temp_drum_*temp_drum_*temp_drum_)
                                                     )
                            / (
                               + (1.0/(4.0*MathConst::MY_PI*rad1i*rad1i*(view_factor_wall_+SMALL_VIEW_FACTOR)))
                              );


						heatFlux[i] -= flux_wall_;				
					}

            fix_raddirectionalHeatFlux->do_forward_comm();
			fix_accumulativeviewfactor->do_forward_comm();
                
    } //end particle loop

	//correct flux because accumulative view factor can be above 1
	for (ii = 0; ii < inum; ii++) 
    {
		i = ilist[ii];
		double procent_over_flux_ = accumulativeviewfactor[i]-1.0;
		if(procent_over_flux_ <= 0.0)
		{
			heatFlux[i] = ((accumulativeviewfactor[i]-1.0)*heatFlux[i])+heatFlux[i];
		}
	}
    fix_heatFlux->do_reverse_comm();

    } //end if-solve_
    solve_counter_ = solve_counter_-1;
}

/* ---------------------------------------------------------------------- */

void FixHeatRad::final_integrate()
{

}

/* ---------------------------------------------------------------------- */

double FixHeatRad::FindDistanceToSegment(double x1, double y1, double z1, double x2, double y2, double z2, double pointX, double pointY, double pointZ)
{
    double diffX = x2 - x1;
    double diffY = y2 - y1;
    double diffZ = z2 - z1;

    if ((diffX == 0) && (diffY == 0) && (diffZ == 0))
    {
        diffX = pointX - x1;
        diffY = pointY - y1;
        diffZ = pointZ - z1;
        return sqrt(diffX * diffX + diffY * diffY +  diffZ * diffZ);
    }

    t = ((pointX - x1) * diffX + (pointY - y1) * diffY + (pointZ - z1) * diffZ) / (diffX * diffX + diffY * diffY + diffZ * diffZ);

    if (t < 0)
    {
        //point is nearest to the first point i.e x1 and y1
        diffX = pointX - x1;
        diffY = pointY - y1;
        diffZ = pointZ - z1;
    }
    else if (t > 1)
    {
        //point is nearest to the end point i.e x2 and y2
        diffX = pointX - x2;
        diffY = pointY - y2;
        diffZ = pointZ - z2;
    }
    else
    {
        //if perpendicular line intersect the line segment.
        diffX = pointX - (x1 + t * diffX);
        diffY = pointY - (y1 + t * diffY);
        diffZ = pointZ - (z1 + t * diffZ);
    }

    //returning shortest distance
    return sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ);
}                             
