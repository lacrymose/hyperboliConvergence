
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <vector>

# include <cassert>

// overload with return value (must be used to construct vector to use RVO)
   template<LawType Law, ImplementedVarSet SolVarT, FluxFunctor<Law> Flux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT>
   std::vector<FluxResult<Law,1,Real>> residualCalc( const std::vector<geom::Volume<1,Real>>& cells,
                                                     const std::vector<geom::Point< 1,Real>>& nodes,
                                                     const Flux flux,
                                                     const Species<Law,Real>& species,
                                                     const SolVarT& ql0,
                                                     const std::vector<SolVarT>& q )
  {
      using FluxRes = FluxResult<Law,1,Real>;
      std::vector<FluxRes> res(q.size());
      residualCalc( cells, nodes, flux, species, ql0, q, res );
      return res;
  }

   template<LawType Law, ImplementedVarSet SolVarT, FluxFunctor<Law> Flux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT>
   void residualCalc( const std::vector<geom::Volume<1,Real>>& cells,
                      const std::vector<geom::Point< 1,Real>>& nodes,
                      const Flux flux,
                      const Species<Law,Real>& species,
                      const SolVarT& ql0,
                      const std::vector<SolVarT>& q,
                            std::vector<FluxResult<Law,1,Real>>& res )
  {
      using FluxRes = FluxResult<Law,1,Real>;
      using Surface = geom::Surface<1,Real>;

   // check mesh sizes match
      const size_t nc = cells.size();

      assert( nc+1 == nodes.size() );
      assert( nc   == res.size() );
      assert( nc   ==   q.size() );

   // reset residual
      for( FluxRes& fr : res ){ fr = FluxRes{}; }

   // accumulate cell residual contributions from the flux across each face
      for( size_t i=0; i<nc-1; i++ )
     {
         const SolVarT ql=q[i];
         const SolVarT qr=q[i+1];

         const Surface face = surface( nodes[i] );

         const FluxRes fr = flux( species, face,  ql, qr );

         res[i]  -=fr;
         res[i+1]+=fr;
     }

   // left boundary : dirichlet condition at initial state
     {
         const size_t i=0;
         const SolVarT qr=q[i];
         const Surface face = surface( nodes[i] );
         const FluxRes fr = flux( species, face, ql0, qr );
         res[i]+=fr;
     }

   // right boundary : outflow, upwind flux
     {
         const size_t i=nc-1;
         const SolVarT ql=q[i];
         const Surface face = surface( nodes[i] );
         const FluxRes fr = flux( species, face, ql, ql );
         res[i]-=fr;
     }
      return;
  }
