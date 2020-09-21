
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

   template<LawType            Law,
            typename       Limiter,
            FluxFunctor<Law>  Flux>
   auto make_muscl_flux( const Limiter& limiter,
                         const Flux&       flux )
  {
      return [limiter,flux]
             <int                      nDim,
              floating_point           Real,
              ImplementedVarSet   SolVarSet,
              ImplementedVarDelta SolVarDel>
            ( const Species<Law,Real>&    species,
              const geom::Surface<nDim,Real> face,
              const SolVarDel&              gradl,
              const SolVarDel&              gradr,
              const SolVarSet&                 ql,
              const SolVarSet&                 qr ) -> FluxResult<Law,nDim,Real>
//       requires ConsistentTypes<Law,nDim,Real,SolVarSet,SolVarDel>
     {
      // central and l/r biased differences in solution basis
         const SolVarDel dqc = qr - ql;
         const SolVarDel dql = gradl - dqc;
         const SolVarDel dqr = gradr - dqc;

      // limited differences in solution basis aligned with background coordinate system
         const SolVarDel slopel = limiter( dqc, dql );
         const SolVarDel sloper = limiter( dqc, dqr );

      // inviscid flux
         return flux( species, face, ql+0.5*slopel,
                                     qr-0.5*sloper );
     };
  }
