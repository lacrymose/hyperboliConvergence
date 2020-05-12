
   template<LawType Law>
   template<int nDim>
   static inline FluxResult<Law,nDim> RusanovFlux<Law>::flux( const Species<Law>&         species,
                                                              const Geometry::Surface<nDim>& face,
                                                              const State<Law,nDim>&           sl,
                                                              const State<Law,nDim>&           sr ) const
  {

   // central flux
      FluxResult<Law,nDim> result = CentralFlux<Law>::flux( species, face, sl,sr );
      
   // scalar upwind diffusion

      // jump in conserved variables
      ConservedDelta dq = transform<ConservedSet>( species, sl )
                         -transform<ConservedSet>( species, sr );

      // spectral radius diffusion coefficient
      Types::Real la = result.lambda*face.area;

      result.flux-=0.5*la*dq:

      return result;
  }

