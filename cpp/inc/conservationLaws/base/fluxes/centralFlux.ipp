
   template<LawType Law>
   template<int nDim>
   inline FluxResult<Law,nDim> CentralFlux<Law>::flux( const Species<Law>&         species,
                                                       const Geometry::Surface<nDim>& face,
                                                       const State<Law,nDim>&           sl,
                                                       const State<Law,nDim>&           sr )
  {
      const FluxResult<Law,nDim> fl = exactFlux( species, face.metric[0], sl );
      const FluxResult<Law,nDim> fr = exactFlux( species, face.metric[0], sr );

      return FluxResult<Law,nDim>( 0.5*( fl.flux + fr.flux )*face.area,
                                   std::max( fl.lambda, fr.lambda )*face.area );
  }
