
# pragma once

# include <spatial/boundary/periodicLoop.h>
# include <spatial/boundary/ghostCellLoop.h>
# include <spatial/boundary/UpdateLoop.h>

# include <solutionField/solutionField.h>

# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/scalarAdvection/boundaryConditions.h>

# include <mesh/mesh.h>

# include <parallalg/array.h>

# include <cassert>

/*
 * Accumulate cell residuals from fluxes over boundary cell faces
 */
   template<int                   nDim,
            ScalarVarSet       SolVarT,
            ScalarVarDelta     SolDelT,
            typename           FluxRes,
            typename     HighOrderFlux,
            floating_point        Real>
      requires   ConsistentTypes<LawType::ScalarAdvection,
                                 nDim,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
   void boundaryResidual( const Mesh<nDim,Real>&                           mesh,
                          const HighOrderFlux&                           hoflux,
                          const Species<LawType::ScalarAdvection,Real>& species,
                          const SolutionField<SolVarT,nDim>&                  q,
                          const par::Array<   SolDelT,nDim>&                 dq,
                                par::Array<   FluxRes,nDim>&                res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );

   // loop over boundaries
      for( unsigned int boundaryId=0; boundaryId<q.nBoundaries; ++boundaryId )
     {
         bool matched=false;
         switch(q.bcTypes[boundaryId])
        {
            case ScalarAdvectionBCs::Periodic :
           {
               matched=true;

               // if (1,3,5), skip (flux already calculated)
               if( boundaryId%2==1 ){ break; }

               // if (0,2,4) check consistent and calculate flux
               assert( q.bctypes[boundaryId+1] == BoundaryType<Law>::Periodic );
               periodicBoundary( mesh, hoflux, boundaryId, q,dq, res );
               break;
           }
            case ScalarAdvectionBCs::RiemannInvariant :
           {
               matched=true;
               ghostBoundary( mesh, hoflux, boundaryId, q,dq, res );
               break;
           }
        }
         assert( matched==true && "boundary has invalid boundary condition flag" );
     }
      return;
  }

/*
 * Update values of boundary cells
 */
   template<int                   nDim,
            ScalarVarSet       SolVarT,
            floating_point        Real>
      requires   ConsistentTypes<LawType::ScalarAdvection,
                                 nDim,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
   void boundaryUpdate( const Mesh<nDim,Real>&                           mesh,
                        const Species<LawType::ScalarAdvection,Real>& species,
                              SolutionField<SolVarT,nDim>&                  q )
  {
   // check array sizes
      assert( mesh.cells.shape() == q.interior.shape() );

      for( unsigned int boundaryId=0; boundaryId<q.nBoundaries; ++boundaryId )
     {
         switch(q.bcTypes[boundaryId])
        {
            case ScalarAdvection::RiemannInvariant :
           {
               using tag = std::integral_constant<ScalarAdvectionBCs,
                                                  ScalarAdvectionBCs::RiemannInvariant>;
               boundaryUpdate( tag{}, mesh, species, boundaryId, q );
               break;
           }
        }
     }

      return;
  }
