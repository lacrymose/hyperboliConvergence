

# pragma once

# include <conservationLaws/base/base.h>

# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <utils/concepts.h>

# include <array>
# include <tuple>

# include <cassert>

/*
 * Accumulate cell residuals from fluxes over all cell faces
 */
   template<LawType                   Law,
            int                      nDim,
            int                    nDimBC,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,nDim,Real,SolVarT,SolDelT>
              && nDimBC==std::max(1,nDim-1)
              && std::is_same_v<FluxRes,fluxresult_t<SolVarT>>
   void residualCalc( const Mesh<nDim,Real>&                        mesh,
                      const HighOrderFlux&                        hoflux,
                      const Species<Law,Real>&                   species,
                      const std::array<BoundaryType<Law>,2*nDim> bctypes,
                      const par::Array<SolVarT,nDim>&                  q,
                      const par::Array<SolDelT,nDim>&                 dq,
                      const par::Array<SolVarT,nDimBC>&               qb,
                            par::Array<FluxRes,nDim>&                res )

  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );

      par::fill( res, FluxRes{} );

   // accumulate cell residual contributions from interior faces
      interiorResidualCalc( mesh, hoflux, species, bctypes, q,dq, res );

   // accumulate cell residual contributions from boundary faces
      boundaryResidualCalc( mesh, hoflux, species, bcfluxes, bctypes, q,dq,qb, res );

      return;
  }


/*
 * Accumulate cell residuals from fluxes over boundary cell faces
 */
   template<LawType                   Law,
            int                      nDim,
            int                    nDimBC,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,nDim,Real,SolVarT,SolDelT>
              && nDimBC==std::max(1,nDim-1)
              && std::is_same_v<FluxRes,fluxresult_t<SolVarT>>
   void boundaryResidualCalc( const Mesh<nDim,Real>&                        mesh,
                              const HighOrderFlux&                        hoflux,
                              const Species<Law,Real>&                   species,
                              const std::array<BoundaryType<Law>,2*nDim> bctypes,
                              const par::Array<SolVarT,nDim>&                  q,
                              const par::Array<SolDelT,nDim>&                 dq,
                              const par::Array<SolVarT,nDimBC>&               qb,
                                    par::Array<FluxRes,nDim>&                res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );

   // loop over boundaries
      for( unsigned int boundaryId=0; boundaryId<2*nDim; boundaryId++ )
     {
         switch(bctypes[boundaryId])
        {
            case BoundaryType<Law>::Periodic :
           {
               if( boundaryId%2==0 )// if (0,2,4) check consistent and calculate flux
              {
                  assert( bctypes[boundaryId+1] == BoundaryType<Law>::Periodic )
                  periodicResidualCalc( mesh, hoflux, species, boundaryId, q,dq,res );
              }
               else{ break; }// if odd, skip (flux already calculated)
           }
            case BoundaryType<Law>::Characteristic :
           {
               characteristicResidualCalc( mesh, species, boundaryId, q,dq,qb, res );
           }
            default : // law specific boundary conditions
           {
               boundaryResidualCalc( mesh, species, bctypes[boundaryId], boundaryId, 
                                     dq,q,qb, res );
           }
        }
     }
      return;
  }
