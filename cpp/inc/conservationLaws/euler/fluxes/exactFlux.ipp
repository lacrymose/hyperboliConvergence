
# include <cmath>

   template<int nDim, floating_point Real>
   FluxResult<LawType::Euler,nDim,Real> exactFlux( const Species<LawType::Euler,Real>&   species,
                                                   const Geometry::Direction<nDim,Real>&  normal,
                                                   const State<LawType::Euler,nDim,Real>&  state )
  {
   // unpack state
      const Real r = state.density();
      const Real p = state.pressure();
      const Real h = state.specificTotalEnthalpy();
      const Real a = std::sqrt( state.speedOfSound2() );

   // face normal velocity and mass flux
      const Real un = projectedVelocity( normal, state );
      const Real mn = r*un;

   // assemble flux
      FluxResult<LawType::Euler,nDim,Real> fr;

   // momentum flux
      for( int i=0; i<nDim; i++ )
     {
         fr.flux[i] = mn*state.velocity(i) + p*normal[i];
     }

   // density and total energy fluxes
      fr.flux[nDim]   = mn;
      fr.flux[nDim+1] = mn*h;

   // max wavespeed
      fr.lambda = std::abs( un ) + a;

      return fr;
  }

