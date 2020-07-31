
# pragma once

# include <conservationLaws/base/base.h>

# include <mesh/mesh.h>

# include <parallalg/array.h>

# include <array>
# include <tuple>

# include <cassert>

   template<LawType                   Law,
            int                      nDim,
            int                    nDimBC,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename...    BoundaryFluxes,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,nDim,Real,SolVarT,SolDelT>
              && nDimBC==std::max(1,nDim-1)
              && std::is_same_v<FluxRes,FluxResult<Law,nDim,Real>>
   void invoke_bc( const HigherOrderFlux&                  hoflux,
                   const std::tuple<BoundaryFluxes...>&  bcfluxes,
                   const size_t                        boundaryId,
                   const BoundaryType<Law>                 bctype,
                   const Mesh<nDim,Real>&                    mesh,
                   const par::Array<SolVarT,nDim>&              q,
                   const par::Array<SolDelT,nDim>&             dq,
                   const par::Array<SolVarT,nDim>&             qb,
                         par::Array<FluxRes,nDim>&            res )
  {
   // check sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );

      qb.get<Indices>()[{i,j}]...

   // accumulate relevant boundary flux on internal residual
      else if( bctype==BoundaryType<Law>::Periodic )
     {
         invoke_periodic_bc( hoflux, boundaryId, mesh, q,dq, res );
     }

   // calculate boundary residual for characteristic boundary conditions
      if( bctype == BoundaryType<Law>::Characteristic )
     {
         invoke_characteristic_bc( hoflux, boundaryId, mesh, q,dq,qb, res );
     }

   // the rest
      ((invoke_functor_bc( bcfluxes.get<Indices>(), Indices, bctype,
                           boundaryId, mesh, q,dq,qb, res )),...);
  }
