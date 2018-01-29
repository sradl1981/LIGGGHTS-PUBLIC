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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2018-     Stefan Radl, Graz University of Technology, Graz
------------------------------------------------------------------------- */

#define COHESION_MODEL_LUBRICATION_VERBOSE true

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_LUBRICATION,lubrication,9)
#else

#ifndef COHESION_MODEL_LUBRICATION_H_
#define COHESION_MODEL_LUBRICATION_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>
#include <algorithm>
#include "math_extra_liggghts.h"
#include "global_properties.h"
#include "fix_property_atom.h"
#include "neighbor.h"

namespace MODEL_PARAMS
{
    inline static ScalarProperty* createMinSeparationDistanceRatioLubrication(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* minSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "minSeparationDistanceRatio", caller);
      return minSeparationDistanceRatioScalar;
    }

    inline static ScalarProperty* createMaxSeparationDistanceRatioLubrication(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* maxSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "maxSeparationDistanceRatio", caller);
      return maxSeparationDistanceRatioScalar;
    }

    inline static ScalarProperty* createModelSwitchLubrication(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* modelSwitchLubrication = MODEL_PARAMS::createScalarProperty(registry, "modelSwitchLubrication", caller);
      return modelSwitchLubrication;
    }

    inline static ScalarProperty* createFluidViscosityLubrication(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidViscosity", caller);
      return fluidViscosityScalar;
    }
}

namespace LIGGGHTS {

namespace ContactModels {

  template<>
  class CohesionModel<COHESION_LUBRICATION> : public CohesionModelBase {

  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *cmb) :
      CohesionModelBase(lmp, hsetup, cmb),
      minSeparationDistanceRatio(0.0),
      maxSeparationDistanceRatio(0.0),
      fluidViscosity(0.),
      modelSwitchLubrication(0),
      modelSwitchHighOrderTermsOn(false), 
      modelSwitchTangentialOn(false),
      history_offset(0)
    {
      history_offset = hsetup->add_history_value("contflag", "0");
      
      if(cmb->is_wall())
        error->warning(FLERR,"Using cohesion model lubrication for walls only supports dry walls");
    }

    void registerSettings(Settings& settings)
    {
        settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("fluidViscosity", &MODEL_PARAMS::createFluidViscosityLubrication);
      registry.registerProperty("minSeparationDistanceRatio", &MODEL_PARAMS::createMinSeparationDistanceRatioLubrication);
      registry.registerProperty("maxSeparationDistanceRatio", &MODEL_PARAMS::createMaxSeparationDistanceRatioLubrication);
      registry.registerProperty("modelSwitchLubrication", &MODEL_PARAMS::createModelSwitchLubrication);

      registry.connect("fluidViscosity", fluidViscosity,"cohesion_model lubrication");
      registry.connect("minSeparationDistanceRatio", minSeparationDistanceRatio,"cohesion_model lubrication");
      registry.connect("maxSeparationDistanceRatio", maxSeparationDistanceRatio,"cohesion_model lubrication");
      double dummyDoubleValue;
      registry.connect("modelSwitchLubrication", dummyDoubleValue, "cohesion_model lubrication");
      modelSwitchLubrication = (int)dummyDoubleValue;
      if(modelSwitchLubrication==1 || modelSwitchLubrication==11) modelSwitchHighOrderTermsOn = true; 
      if(modelSwitchLubrication>9) modelSwitchTangentialOn = true;
      ln1overMinSeparationDistanceRatio = log(1./minSeparationDistanceRatio);

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model lubrication cannot work with coarsegraining at this stage");

      // error check on contact distance factor
      neighbor->register_contact_dist_factor(maxSeparationDistanceRatio); 
      if(maxSeparationDistanceRatio < 1.1)
            error->one(FLERR,"\n\ncohesion model lubrication requires maxSeparationDistanceRatio >= 1.1. Please increase this value, otherwise the lubrication model may result in incorrect results.\n");
    }

    inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const double radi = sidata.radi;
      const double radj = sidata.is_wall ? radi : sidata.radj;

      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;
      double * const contflag = &sidata.contact_history[history_offset];
      // store for noCollision
      contflag[0] = 1.0;

      const double rEff = radi*radj / (radi+radj); //for equal spheres this is rad/2
    
      // viscous normal force
      const double stokesPreFactor = -6.*M_PI*fluidViscosity*rEff;
      const double FviscN = modelSwitchHighOrderTermsOn ?
                            stokesPreFactor
                          * sidata.vn
                          * (sidata.is_wall ?
                                particleWallNormalHighOrder(minSeparationDistanceRatio) 
                               :particleParticleNormalHighOrder(minSeparationDistanceRatio) )
                          : 
                            stokesPreFactor
                          * sidata.vn
                          * (sidata.is_wall ?
                                particleWallNormalSimple(minSeparationDistanceRatio) 
                               :particleParticleNormalSimple(minSeparationDistanceRatio));


      // viscous tangential force 
      const double FviscT_over_vt = modelSwitchTangentialOn ?
                                    stokesPreFactor
                                  * particleParticleTangentialVT(minSeparationDistanceRatio)
                                  :          
                                    0;

      const double Ft1 = FviscT_over_vt    * sidata.vtr1; //Simplified treatment since only total sliding speed is known!
      const double Ft2 = FviscT_over_vt    * sidata.vtr2;
      const double Ft3 = FviscT_over_vt    * sidata.vtr3;

      // torques
      const double tor1 = sidata.en[1] * Ft3 - sidata.en[2] * Ft2;
      const double tor2 = sidata.en[2] * Ft1 - sidata.en[0] * Ft3;
      const double tor3 = sidata.en[0] * Ft2 - sidata.en[1] * Ft1;

      // add to fn, Ft
      if(tangentialReduce_) sidata.Fn += FviscN;  

      // apply normal and tangential force
      const double fx = FviscN * sidata.en[0] + Ft1;
      const double fy = FviscN * sidata.en[1] + Ft2;
      const double fz = FviscN * sidata.en[2] + Ft3;

#if COHESION_MODEL_LUBRICATION_VERBOSE
      printf("Coh_lubrication_surfacesIntersect, FviscN: %.3g, vn: %.3g \n", 
                FviscN, sidata.vn
            );
#endif

      // return resulting forces
      if(sidata.is_wall) {
        const double area_ratio = sidata.area_ratio;
        i_forces.delta_F[0] += fx * area_ratio;
        i_forces.delta_F[1] += fy * area_ratio;
        i_forces.delta_F[2] += fz * area_ratio;
        i_forces.delta_torque[0] += -sidata.cri * tor1 * area_ratio;
        i_forces.delta_torque[1] += -sidata.cri * tor2 * area_ratio;
        i_forces.delta_torque[2] += -sidata.cri * tor3 * area_ratio;
      } else {
        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;
        i_forces.delta_torque[0] += -sidata.cri * tor1;
        i_forces.delta_torque[1] += -sidata.cri * tor2;
        i_forces.delta_torque[2] += -sidata.cri * tor3;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;
        j_forces.delta_torque[0] += -sidata.crj * tor1;
        j_forces.delta_torque[1] += -sidata.crj * tor2;
        j_forces.delta_torque[2] += -sidata.crj * tor3;
      }
    }

    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {

      double * const contflag = &scdata.contact_history[history_offset];

      const int i = scdata.i;
      const int j = scdata.j;
      const double radi = scdata.radi;
      const double radj = scdata.is_wall ? radi : scdata.radj;
      const double r = sqrt(scdata.rsq);

#if COHESION_MODEL_LUBRICATION_VERBOSE
      printf("Coh_lubrication_surfacesClose, radi/radj: %.3g / %.3g, r: %.3g, critical distance: %.3g \n", 
                radi, radj, r, maxSeparationDistanceRatio*(radi+radj)
            );
#endif

      //exit in case particles are too far apart
      if(r>maxSeparationDistanceRatio*(radi+radj))
      {
        contflag[0] = 0.0;
        scdata.has_force_update = false;
        return;
      }

      const double dist = scdata.is_wall ? r - radi : r - (radi + radj);
      const double rEff = radi*radj / (radi+radj);

      double **v = atom->v;

      //reset contact flag (spare)
      contflag[0] = 0.0;

      // calculate vn and vt since not in struct
      const double rinv = 1.0 / r;
      const double dx = scdata.delta[0];
      const double dy = scdata.delta[1];
      const double dz = scdata.delta[2];
      const double enx = dx * rinv; //unity vector
      const double eny = dy * rinv;
      const double enz = dz * rinv;

      // relative translational velocity
      const double vr1 = v[i][0] - v[j][0];
      const double vr2 = v[i][1] - v[j][1];
      const double vr3 = v[i][2] - v[j][2];

      // normal component
      const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
      const double vn1 = vn * enx;
      const double vn2 = vn * eny;
      const double vn3 = vn * enz;

      // tangential component (translational motion ONLY!)
      const double vt1 = vr1 - vn1;
      const double vt2 = vr2 - vn2;
      const double vt3 = vr3 - vn3;

      // relative rotational velocity
      double wr1, wr2, wr3;
      double const *omega_i = atom->omega[i];
      double const *omega_j = atom->omega[j];

      if(scdata.is_wall) {
            wr1 = radi * omega_i[0] * rinv;
            wr2 = radi * omega_i[1] * rinv;
            wr3 = radi * omega_i[2] * rinv;
      } else {
            wr1 = radi * omega_i[0] + radj * omega_j[0];
            wr2 = radi * omega_i[1] + radj * omega_j[1];
            wr3 = radi * omega_i[2] + radj * omega_j[2];
      }

      // relative velocities (Rotational motion ONLY!)
      const double vtRot1 = -(enz * wr2 - eny * wr3);
      const double vtRot2 = -(enx * wr3 - enz * wr1);
      const double vtRot3 = -(eny * wr1 - enx * wr2);

      // viscous normal force
      const double stokesPreFactor = -6.*M_PI*fluidViscosity*rEff;
      const double FviscN = modelSwitchHighOrderTermsOn ?
                            stokesPreFactor
                          * vn
                          * (scdata.is_wall ?
                                particleWallNormalHighOrder(std::max(minSeparationDistanceRatio,dist/rEff)) 
                               :particleParticleNormalHighOrder(std::max(minSeparationDistanceRatio,dist/rEff)) )
                          : 
                            stokesPreFactor
                          * vn
                          * (scdata.is_wall ?
                                particleWallNormalSimple(std::max(minSeparationDistanceRatio,dist/rEff)) 
                               :particleParticleNormalSimple(std::max(minSeparationDistanceRatio,dist/rEff)));


      // viscous tangential force 
      const double FviscT_over_vt = modelSwitchTangentialOn ?
                                    stokesPreFactor
                                  * (scdata.is_wall ?
                                        particleWallTangentialVT(std::max(minSeparationDistanceRatio,dist/rEff))
                                       :particleParticleTangentialVT(std::max(minSeparationDistanceRatio,dist/rEff)))
                                  :          
                                    0;

      const double FviscT_over_vtRot = modelSwitchTangentialOn ?
                                       stokesPreFactor
                                     * (scdata.is_wall ?
                                        particleWallTangentialVTRot(std::max(minSeparationDistanceRatio,dist/rEff))
                                       :particleParticleTangentialVTRot(std::max(minSeparationDistanceRatio,dist/rEff)))
                                     :          
                                       0;

      const double Ft1 = FviscT_over_vt    * vt1
                       + FviscT_over_vtRot * vtRot1;
      const double Ft2 = FviscT_over_vt    * vt2
                       + FviscT_over_vtRot * vtRot2;
      const double Ft3 = FviscT_over_vt    * vt3
                       + FviscT_over_vtRot * vtRot3;

      // torques (EXCLUDES twisting torque1)
      const double tor1 = eny * Ft3 - enz * Ft2;
      const double tor2 = enz * Ft1 - enx * Ft3;
      const double tor3 = enx * Ft2 - eny * Ft1;

      // apply normal and TOTAL tangential force
      const double fx = FviscN * enx + Ft1;
      const double fy = FviscN * eny + Ft2;
      const double fz = FviscN * enz + Ft3;

      scdata.has_force_update = true;

#if COHESION_MODEL_LUBRICATION_VERBOSE
      printf("Coh_lubrication_surfacesClose, FviscN: %.3g, vn: %.3g, FviscT: %.3g %.3g %.3g \n", 
                FviscN, vn,
                Ft1, Ft2, Ft3
            );
#endif


      // return resulting forces
      if(scdata.is_wall) {
            const double area_ratio = scdata.area_ratio;
            i_forces.delta_F[0] += fx * area_ratio;
            i_forces.delta_F[1] += fy * area_ratio;
            i_forces.delta_F[2] += fz * area_ratio;
            i_forces.delta_torque[0] += -radi * tor1 * area_ratio;
            i_forces.delta_torque[1] += -radi * tor2 * area_ratio;
            i_forces.delta_torque[2] += -radi * tor3 * area_ratio;
      } else {
            i_forces.delta_F[0] += fx;
            i_forces.delta_F[1] += fy;
            i_forces.delta_F[2] += fz;
            i_forces.delta_torque[0] += -radi * tor1; // using radius here, not contact radius
            i_forces.delta_torque[1] += -radi * tor2;
            i_forces.delta_torque[2] += -radi * tor3;

            j_forces.delta_F[0] -= fx;
            j_forces.delta_F[1] -= fy;
            j_forces.delta_F[2] -= fz;
            j_forces.delta_torque[0] += -radj * tor1; // using radius here, not contact radius
            j_forces.delta_torque[1] += -radj * tor2;
            j_forces.delta_torque[2] += -radj * tor3;
      }
    }

  private:
    double minSeparationDistanceRatio, maxSeparationDistanceRatio, fluidViscosity;
    int    modelSwitchLubrication;
    bool   modelSwitchHighOrderTermsOn, modelSwitchTangentialOn;
    double ln1overMinSeparationDistanceRatio;
    int history_offset;
    bool tangentialReduce_;

    //-------------------------------------------------------------------
    //Functions implementing the force coefficient calculations PARTICLE-PARTICLE
    //NORMAL FORCE
    double particleParticleNormalSimple(double normalizedSepDistance)    
    {
        //Dance and Maxey, J Comp Phys 189:212-238  & Kosek et al., Langmuir 32, p.8453. Only first term! 
        //Note, we use rEff for the Stokes prefactor AND normalizedSepDistance: so need to multiply with 4! (see also Brenner, 1961, Chem Eng Sci
        return 1.0/normalizedSepDistance;
    }
    double particleParticleNormalHighOrder(double normalizedSepDistance)    
    {
        //Dance and Maxey, J Comp Phys 189:212-238  & Kosek et al., Langmuir 32, p.8453. Full terms. Note, we use rEff, so need to multiply with 4 and 2 
        return  1.0  / normalizedSepDistance 
              - 0.45 * log(0.5*normalizedSepDistance)
              - 2.6786e-2 * normalizedSepDistance * log(0.5*normalizedSepDistance);
    }

    //TANGENIAL FORCE
    double particleParticleTangentialVT(double normalizedSepDistance)    
    {
        //Dance and Maxey, J Comp Phys 189:212-238 & Kosek et al., Langmuir 32, p.8453. Only first term! 
        //Note, we use rEff for the Stokes prefactor AND normalizedSepDistance: so need to multiply with 2! 
        return -0.333333333333 * log(0.5*normalizedSepDistance);
    }
    double particleParticleTangentialVTRot(double normalizedSepDistance)    
    {
        //Dance and Maxey, J Comp Phys 189:212-238  & Kosek et al., Langmuir 32, p.8453. Only first term! 
        //Note, we use rEff for the Stokes prefactor AND normalizedSepDistance: so need to multiply with 2 and 1!
        return -0.333333333333 * log(0.5*normalizedSepDistance)
               -0.083333333333 * normalizedSepDistance * log(0.5*normalizedSepDistance);
    }


    //-------------------------------------------------------------------
    //Functions implementing the force coefficient calculations PARTICLE-WALL
    double particleWallNormalSimple(double normalizedSepDistance)    
    {
        //Dance and Maxey, J Comp Phys 189:212-238 & Kosek et al., Langmuir 32, p.8453. First term only
        return  1.0  / normalizedSepDistance;
    }
    double particleWallNormalHighOrder(double normalizedSepDistance)    
    {
        //Dance and Maxey, J Comp Phys 189:212-238 & Kosek et al., Langmuir 32, p.8453. Full terms.
        return  1.0  / normalizedSepDistance 
              - 0.20 * log(normalizedSepDistance)
              - 4.7619e-2 * normalizedSepDistance * log(normalizedSepDistance);
    }

    //TANGENIAL FORCE
    double particleWallTangentialVT(double normalizedSepDistance)    
    {
        //Dance and Maxey, J Comp Phys 189:212-238 & Kosek et al., Langmuir 32, p.8453. Full terms.
        return -0.533333333333 * log(normalizedSepDistance)
               -0.170666667 * normalizedSepDistance* log(normalizedSepDistance);
    }
    double particleWallTangentialVTRot(double normalizedSepDistance)    
    {
        //Dance and Maxey, J Comp Phys 189:212-238 & Kosek et al., Langmuir 32, p.8453. Full terms.
        return -0.133333333333 * log(normalizedSepDistance)
               -0.229333333333 * normalizedSepDistance* log(normalizedSepDistance);
    }

  };
}
}
#endif // COHESION_MODEL_CAPILLARY_H_
#endif
