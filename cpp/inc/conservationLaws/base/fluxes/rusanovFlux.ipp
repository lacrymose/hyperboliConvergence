
   template<LawType Law>
   template<int nDim>
   inline FluxResult<Law,nDim> RusanovFlux<Law>::flux( const Species<Law>&         species,
                                                       const Geometry::Surface<nDim>& face,
                                                       const State<Law,nDim>&           sl,
                                                       const State<Law,nDim>&           sr )
  {

   // types needed for upwind diffusion calculation
      using ConservedSet   = VariableSet<  Law,nDim,BasisType<Law>::Conserved>;
      using ConservedDelta = VariableDelta<Law,nDim,BasisType<Law>::Conserved>;

   // central flux
      FluxResult<Law,nDim> result = CentralFlux<Law>::flux( species, face, sl,sr );
      
   // scalar upwind diffusion

      // jump in conserved variables
      const ConservedDelta dq = state2Set<ConservedSet>( species, sr )
                               -state2Set<ConservedSet>( species, sl );

      // spectral radius diffusion coefficient
      const Types::Real la = result.lambda*face.area;

      result.flux-=0.5*la*dq;

      return result;
  }

