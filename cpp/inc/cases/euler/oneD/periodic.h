
# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <parallalg/array.h>

# include <cmath>

   enum struct Periodic1D
  {
      Soundwave,
      AcousticEntropy
  };

//   template<EulerVarSet VarT, floating_point Real>
//      requires SameFPType<VarT,Real>
//   par::Array<VarT,1> periodic_initial_solution( const Periodic1D problem,
//                                                 const Species<LawType::Euler,Real>& species,
//                                                 const par::Array<geom::Volume<1,Real>,1>& cells )
//  {
//      switch( problem )
//     {
//         case( Periodic1D::Soundwave ) :
//        {
//            return soundwave_initial_solution<VarT>( species, cells );
//        }
//     }
//  }

/*
 * simple forward-travelling sinusoidal soundwave
 *    nwavelength periods over a domain of length
 *    u0,rho0 = 1, p0 = 1./mach^2
 */
     template<EulerVarSet VarT, floating_point Real>
        requires SameFPType<VarT,Real>
     par::Array<VarT,1> soundwave_initial_solution( const Species<LawType::Euler,Real>&     species,
                                                    const par::Array<geom::Volume<1,Real>,1>& cells,
                                                    const Real               length,
                                                    const Real                 mach,
                                                    const Real            amplitude,
                                                    const unsigned int nwavelengths )
  {
      using PrimitiveVarT = VariableSet<LawType::Euler,1,EulerBases::Primitive,Real>;

      par::Array<VarT,1> q(cells.shape());

   // background state
      const Real u0 = 1.;
      const Real r0 = 1.;
      const Real a0 = u0/mach;
      const Real p0 = r0*a0*a0/species.gamma;

   // wavenumber
      const Real kappa = 2.*M_PI*nwavelengths/length;

      const size_t nc = cells.shape(0);

      for( size_t i=0; i<nc; i++ )
     {
         const Real x = cells[{i}].centre[0];

      // set pressure
         const Real p = p0*( 1. + mach*amplitude*sin( kappa*x ) );

      // calculate density from isentropic relation
         const Real r = r0*pow( p/p0, 1./species.gamma );

      // calculate local speed of sound a^2 = gamma*p/rho
         const Real a = sqrt( species.gamma*p/r );

      // calculate velocity from riemann invariant u-2*a/(gamma-1) for u+a wave
         const Real u = u0 + 2.*( a - a0 )/(species.gamma-1.);

         q[{i}] = set2Set<VarT>( species, PrimitiveVarT{u,r,p} );
     }
      return q;
  }



/*
 * forward-travelling gaussian soundwave travels through gaussian entropy bump
 *    u0,rho0 = 1, p0 = 1./mach^2
 */
     template<EulerVarSet VarT, floating_point Real>
        requires SameFPType<VarT,Real>
     par::Array<VarT,1> acoustic_entropy_initial_solution( const Species<LawType::Euler,Real>&     species,
                                                           const par::Array<geom::Volume<1,Real>,1>& cells,
                                                           const Real               length,
                                                           const Real                 mach,
                                                           const Real             centre_a,
                                                           const Real             centre_s,
                                                           const Real              width_a,
                                                           const Real              width_s,
                                                           const Real          amplitude_a,
                                                           const Real          amplitude_s )
  {
      using PrimitiveVarT = VariableSet<LawType::Euler,1,EulerBases::Primitive,Real>;

      par::Array<VarT,1> q(cells.shape());

   // background state
      const Real u0 = 1.;
      const Real r0 = 1.;
      const Real a0 = u0/mach;
      const Real p0 = r0*a0*a0/species.gamma;

   // gaussian parameters

      // centre
      const Real mu_a = centre_a;
      const Real mu_s = centre_s;

      // width
      const Real sigma_a = width_a;
      const Real sigma_s = width_s;

      // coefficients
      const Real c0_a = amplitude_a / ( sigma_a * sqrt( 2.*M_PI ) );
      const Real c0_s = amplitude_s / ( sigma_s * sqrt( 2.*M_PI ) );

      const size_t nc = cells.shape(0);

      for( size_t i=0; i<nc; i++ )
     {
      // relative location along domain
         const Real x = cells[{i}].centre[0]/length;

      // acoustic and entropy gaussians
         // exponents
         const Real e_a = -0.5*( x-mu_a )*( x-mu_a )/( sigma_a*sigma_a );
         const Real e_s = -0.5*( x-mu_s )*( x-mu_s )/( sigma_s*sigma_s );

         const Real g_a = c0_a * exp( e_a );
         const Real g_s = c0_s * exp( e_s );

      // set pressure
         const Real p = p0*( 1. + mach*g_a );

      // calculate velocity from isentropic acoustic relations
         // isentropic density
         const Real r_i = r0*pow( p/p0, 1./species.gamma );

         // isentropic local speed of sound
         const Real a_i = sqrt( species.gamma*p/r_i );

      // calculate actual density from isentropic relation + entropy displacement
         const Real r = r0*pow( p/p0, 1./species.gamma ) + g_s*mach;

         // local speed of sound
         [[maybe_unused]] const Real a = sqrt( species.gamma*p/r );

      // calculate velocity from riemann invariant u-2*a/(gamma-1) for u+a wave
         const Real u = u0 + 2.*( a_i - a0 )/(species.gamma-1.);

         q[{i}] = set2Set<VarT>( species, PrimitiveVarT{u,r,p} );
     }

      return q;
  }


