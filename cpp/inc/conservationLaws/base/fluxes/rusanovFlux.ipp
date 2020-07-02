
   template<LawType Law>
   template<int nDim, floating_point Real>
   FluxResult<Law,nDim,Real> RusanovFlux<Law>::flux( const Species<Law,Real>&     species,
                                                     const geom::Surface<nDim,Real>& face,
                                                     const State<Law,nDim,Real>&       sl,
                                                     const State<Law,nDim,Real>&       sr )
  {

   // types needed for upwind diffusion calculation
      using ConservedSet   = VariableSet<  Law,nDim,BasisType<Law>::Conserved,Real>;
      using ConservedDelta = VariableDelta<Law,nDim,BasisType<Law>::Conserved,Real>;
      using FluxRes = FluxResult<Law,nDim,Real>;

   // central flux
      FluxRes result = CentralFlux<Law>::flux( species, face, sl,sr );

   // scalar upwind diffusion

      // jump in conserved variables
      const ConservedDelta dq = state2Set<ConservedSet>( species, sr )
                               -state2Set<ConservedSet>( species, sl );

      // spectral radius diffusion coefficient
      const Real la = result.lambda*face.area;

      result.flux-=0.5*la*dq;

      return result;
  }


   template<LawType Law>
   template<int nDim, floating_point Real>
   FluxResult<Law,nDim,Real> RusanovDissipation<Law>::flux( const Species<Law,Real>&     species,
                                                            const geom::Surface<nDim,Real>& face,
                                                            const State<Law,nDim,Real>&       sl,
                                                            const State<Law,nDim,Real>&       sr )
  {

   // types needed for upwind diffusion calculation
      using ConservedSet   = VariableSet<  Law,nDim,BasisType<Law>::Conserved,Real>;
      using ConservedDelta = VariableDelta<Law,nDim,BasisType<Law>::Conserved,Real>;
      using FluxRes = FluxResult<Law,nDim,Real>;

   // spectral radius
      const Real lambda = std::max( spectralRadius( face.metric[0], sl ),
                                    spectralRadius( face.metric[0], sr ) );

   // scalar upwind diffusion

      // jump in conserved variables
      const ConservedDelta dq = state2Set<ConservedSet>( species, sr )
                               -state2Set<ConservedSet>( species, sl );

      // spectral radius diffusion coefficient
      const Real la = lambda*face.area;

      return FluxRes( -0.5*la*dq, lambda );
  }

