
   template<LawType Law>
   template<int nDim>
   static inline FluxResult<Law,nDim> CentralFlux<Law>::flux( const Species<Law>&         species,
                                                              const Geometry::Surface<nDim>& face,
                                                              const State<Law,nDim>&           sl,
                                                              const State<Law,nDim>&           sr ) const
  {
      fl = exactFlux( species, face, sl );
      fr = exactFlux( species, face, sr );
      return FluxResult<Law,nDim>( 0.5*( fl.flux + fr.flux )*face.area,
                                   std::max( fl.lambda, fr.lambda ) );
  }
