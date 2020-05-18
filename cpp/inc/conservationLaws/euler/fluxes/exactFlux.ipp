
# include <cmath>

   template<int nDim>
   FluxResult<LawType::Euler,nDim> exactFlux( const Species<LawType::Euler>&   species,
                                              const Geometry::Direction<nDim>&  normal,
                                              const State<LawType::Euler,nDim>&  state )
  {
   // unpack state
      const Types::Real r = state.density();
      const Types::Real p = state.pressure();
      const Types::Real h = state.specificTotalEnthalpy();
      const Types::Real a = std::sqrt( state.speedOfSound2() );

   // face normal velocity and mass flux
      const Types::Real un = projectedVelocity( normal, state );
      const Types::Real mn = r*un;

   // assemble flux
      FluxResult<LawType::Euler,nDim> fr;

   // momentum flux
      for( int i=0; i<nDim; i++ )
     {
         fr.flux[i] = mn*state.velocity(i) + p*normal[i];
     }

   // density and total energy fluxes
      fr.flux[nDim]   = mn;
      fr.flux[nDim+1] = mn*h;

   // max wavespeed
      fr.lambda= std::abs( un ) + a;

      return fr;
  }

