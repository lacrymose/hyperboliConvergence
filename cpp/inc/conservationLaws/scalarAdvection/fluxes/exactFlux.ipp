
/*
 * exact flux according of the Scalar Advection equation dq/dt + d(uq)/dx = 0
 *    assumes advecting velocity is constant in time (but not necessarily in space), ie:
 *       flux[0:nDim-1] is zero (the advecting velocity) 
 *       flux[nDim] is the projection of uq onto the flux direction
 */
   template<int nDim>
   FluxResult<LawType::ScalarAdvection,nDim> exactFlux( const Species<LawType::ScalarAdvection>&  species,
                                                        const Geometry::Direction<nDim>&           normal,
                                                        const State<LawType::ScalarAdvection,nDim>& state )
  {
      const Types::Real un = projectedVelocity( normal, state );
      const Types::Real q = state.scalar();

      FluxResult<LawType::ScalarAdvection,nDim> result;
      result.flux[nDim]=un*q;
      result.lambda=un;

      return result;
  }

