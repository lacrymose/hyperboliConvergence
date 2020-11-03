
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <limiters/limiter.h>

   template<LawType            Law,
            typename       Limiter,
            FluxFunctor<Law>  Flux>
   auto make_muscl_flux( const Limiter& limiter,
                         const Flux&       flux )
  {
      return [limiter,flux]
             <int                      nDim,
              floating_point           Real,
              ImplementedVarSet   SolVarSet>
            ( const Species<Law,Real>&       species,
              const geom::Surface<nDim,Real>&   face,
              const geom::Volume<nDim,Real>&  cell_l,
              const geom::Volume<nDim,Real>&  cell_r,
              const SolVarSet&                   q_l,
              const SolVarSet&                   q_r,
              const lsq::XMetric<nDim,Real>&   dxm_l,
              const lsq::XMetric<nDim,Real>&   dxm_r,
              const lsq::QMetric<SolVarSet>&   dqm_l,
              const lsq::QMetric<SolVarSet>&   dqm_r ) -> fluxresult_t<SolVarSet>
//       requires ConsistentTypes<Law,nDim,Real,SolVarSet>
     {
         using Direction = geom::Direction<nDim,Real>;
         using SolVarDel = vardelta_t<SolVarSet>;

      // central deltas
         const Direction dx_c = cell_r.centre - cell_l.centre 
         const SolVarDel dq_c = q_r - q_l;

      // distances
         const Real dl = geom::length(   face.centre - cell_l.centre );
         const Real dr = geom::length( cell_r.centre -   face.centre );

      // central gradient
         const SolVarDel dqdx_c = dq_c/(dl+dr);

      // biased left/right gradients
         const std::pair dqdx_lr = lsq::bias_solves( dx_c, dq_c,
                                                     dxm_l,dxm_r,
                                                     dqm_l,dqm_r,
                                                     face.normal );

      // limited left/right gradients
         const SolVarDel slope_l = limiter( dqdx_c, dqdx_lr.first  );
         const SolVarDel slope_r = limiter( dqdx_c, dqdx_lr.second );

      // inviscid flux
         return flux( species, face, ql + dl*slope_l,
                                     qr - dr*slope_r );
     };
  }

   template<LawType            Law,
            FluxFunctor<Law>  Flux>
   auto make_muscl_flux( const Limiters::NoLimit1&,
                         const Flux&          flux )
  {
      return [flux]
             <int                      nDim,
              floating_point           Real,
              ImplementedVarSet   SolVarSet,
              ImplementedVarDelta SolVarDel>
            ( const Species<Law,Real>&     species,
              const geom::Surface<nDim,Real>& face,
              const geom::Volume<nDim,Real>&,
              const geom::Volume<nDim,Real>&,
              const SolVarSet&                 ql,
              const SolVarSet&                 qr,
              const lsq::XMetric<nDim,Real>&,
              const lsq::XMetric<nDim,Real>&,
              const lsq::QMetric<SolVarSet>&,
              const lsq::QMetric<SolVarSet>& ) -> fluxresult_t<SolVarSet>
//       requires ConsistentTypes<Law,nDim,Real,SolVarSet,SolVarDel>
     {
      // first order inviscid flux
         return flux( species, face, ql, qr );
     };
  }

