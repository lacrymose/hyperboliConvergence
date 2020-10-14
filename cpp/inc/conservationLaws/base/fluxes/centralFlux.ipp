
   template<LawType Law>
   template<int nDim, floating_point Real>
   FluxResult<Law,nDim,Real> CentralFlux<Law>::flux( const Species<Law,Real>&     species,
                                                     const geom::Surface<nDim,Real>& face,
                                                     const State<Law,nDim,Real>&       sl,
                                                     const State<Law,nDim,Real>&       sr )
  {
      const FluxResult<Law,nDim,Real> fl = exactFlux( species, face.metric[0], sl );
      const FluxResult<Law,nDim,Real> fr = exactFlux( species, face.metric[0], sr );

      return FluxResult<Law,nDim,Real>( 0.5*( fl.flux + fr.flux )*face.area,
                                        std::max( fl.lambda, fr.lambda )*face.area );
  }
