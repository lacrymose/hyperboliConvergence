
/*
 * exact flux according of the Scalar Advection equation dq/dt + d(uq)/dx = 0
 *    assumes advecting velocity is constant in time (but not necessarily in space), ie:
 *       flux[0:nDim-1] is zero (the advecting velocity) 
 *       flux[nDim] is the projection of uq onto the flux direction
 */
   template<int nDim, floating_point Real>
   FluxResult<LawType::ScalarAdvection,nDim,Real> exactFlux( const Species<LawType::ScalarAdvection,Real>&       species,
                                                             const Geometry::Direction<nDim,Real>&                normal,
                                                             const State<LawType::ScalarAdvection,nDim,Real>&      state )
  {
      const Real un = projectedVelocity( normal, state );
      const Real q  = state.scalar();

      FluxResult<LawType::ScalarAdvection,nDim,Real> result{};
      result.flux[nDim]=un*q;
      result.lambda=un;

      return result;
  }

