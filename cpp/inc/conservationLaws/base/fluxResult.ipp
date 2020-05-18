
/*
 * in-place accumulation for FluxResults
 *    assumes accumulation is taking place for fluxes into a volume:
 *    this means that spectral radii are always added together, even if fluxes are subtracted
 *    in-place accumulation should not be used within interface flux functions for eg adding together central and diffusive fluxes. In this case, the spectral radius is the maximum of the two flux components
 */
   template<LawType Law, int nDim>
   FluxResult<Law,nDim>& FluxResult<Law,nDim>::operator+=( const FluxResult& fr )
  {
      flux   += fr.flux;
      lambda += fr.lambda;
      return *this;
  }

   template<LawType Law, int nDim>
   FluxResult<Law,nDim>& FluxResult<Law,nDim>::operator-=( const FluxResult& fr )
  {
      flux   -= fr.flux;
      lambda += fr.lambda;
      return *this;
  }

