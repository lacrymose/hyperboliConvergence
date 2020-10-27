
# pragma once

# include <spatial/boundary/boundaryCondition.h>

# include <conservationLaws/base/base.h>
# include <solutionField/solutionField.h>

# include <geometry/geometry.h>
# include <mesh/mesh.h>

# include <parallalg/array.h>

# include <cassert>


/*
 * Accumulate cell residual contributions over boundary flux boundary faces for 1D domain
 */
   template<LawType                   Law,
            BoundaryType<Law>      BCType,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename           UpdateFunc,
            typename           BCFluxFunc,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 1,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
   void boundaryResidual( const BoundaryCondition<Law,BCType,UpdateFunc,BCFluxFunc>& bc,
                          const size_t                                       boundaryId,
                          const Mesh<1,Real>&                                      mesh,
                          const HighOrderFlux&                                   hoflux,
                          const Species<Law,Real>&                              species,
                          const SolutionField<SolVarT,1>&                             q,
                          const par::DualArray<std::array<SolDelT,1>,1>&             dq,
                                par::DualArray<FluxRes,1>&                          res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

   // boundary condition type selection was correct?
      assert( q.bcTypes[boundaryId] == BCType );

      auto& boundaryFlux = std::get<1>(bc.funcs);

      using CellIdx = typename SolutionField<SolVarT,1>::VarField::IdxType;
      using NodeIdx = typename Mesh<1,Real>::NodeArray::IdxType;

   // first order boundaries
      if( boundaryId==0 ) // left boundary
     {
         const NodeIdx ip{0};
         const CellIdx ic{0};
         const CellIdx ib0{0};
         const CellIdx ib1{1};

         const FluxRes fr = boundaryFlux( species,
                                          surface( mesh.nodes[ip] ),
                                          mesh.cells[ic],
                                          dq[ic][0],
                                          q.interior[ic],
                                          q.boundary[0][ib0],
                                          q.boundary[0][ib1] );
         res[ic]+=fr;
     }
      else if( boundaryId==1 ) // right boundary
     {
         const size_t nc = q.interior.shape(0);

         const NodeIdx ip{nc-1};
         const CellIdx ic{nc-1};
         const CellIdx ib0{0};
         const CellIdx ib1{1};

         const FluxRes fr = boundaryFlux( species,
                                          flip( surface( mesh.nodes[ip] ) ),
                                          mesh.cells[ic],
                                          dq[ic][0],
                                          q.interior[ic],
                                          q.boundary[0][ib0],
                                          q.boundary[0][ib1] );
         res[ic]+=fr;
     }
      else{ assert( false && "invalid boundary id for 1D ghost cell flux, must be 0 or 1" ); }

      return;
  }

/*
 * Accumulate cell residual contributions over boundary flux boundary faces for 2D domain
 */
   template<LawType                   Law,
            BoundaryType<Law>      BCType,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename           UpdateFunc,
            typename           BCFluxFunc,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 2,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
   void boundaryResidual( const BoundaryCondition<Law,BCType,UpdateFunc,BCFluxFunc>& bc,
                          const size_t                           boundaryId,
                          const Mesh<2,Real>&                          mesh,
                          const HighOrderFlux&                       hoflux,
                          const Species<Law,Real>&                  species,
                          const SolutionField<SolVarT,2>&                 q,
                          const par::DualArray<std::array<SolDelT,2>,2>& dq,
                                par::DualArray<FluxRes,2>&              res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

   // boundary condition type selection was correct?
      assert( q.bcTypes[boundaryId] == BCType );

      const size_t ni = q.interior.shape(0);
      const size_t nj = q.interior.shape(1);

      using CellIdx = typename SolutionField<SolVarT,2>::VarField::IdxType;
      using NodeIdx = typename Mesh<2,Real>::NodeArray::IdxType;

      auto& boundaryFlux = std::get<1>(bc.funcs);

      if( boundaryId==0 ) // left face
     {
         const size_t i=0;
         for( size_t j=0; j<nj; ++j )
        {
         // interior cell  index
            const CellIdx ic{i,j};

         // boundary cell index
            const CellIdx ib0{j,0};
            const CellIdx ib1{j,1};
   
         // face node indices
            const NodeIdx ip0{i,j  };
            const NodeIdx ip1{i,j+1};
   
            const FluxRes fr = boundaryFlux( species,
                                             surface( mesh.nodes[ip0],
                                                      mesh.nodes[ip1] ),
                                             mesh.cells[ic],
                                             dq[ic][0],
                                             q.interior[ic],
                                             q.boundary[0][ib0],
                                             q.boundary[0][ib1] );
            res[ic]+=fr;
        }
     }
      else if( boundaryId==1 ) // right face
     {
         const size_t i=ni-1;
         for( size_t j=0; j<nj; ++j )
        {
         // interior cell  index
            const CellIdx ic{i,j};

         // boundary cell index
            const CellIdx ib0{j,0};
            const CellIdx ib1{j,1};
   
         // face node indices
            const NodeIdx ip0{i+1,j+1};
            const NodeIdx ip1{i+1,j  };
   
            const FluxRes fr = boundaryFlux( species,
                                             surface( mesh.nodes[ip0],
                                                      mesh.nodes[ip1] ),
                                             mesh.cells[ic],
                                             dq[ic][0],
                                             q.interior[ic],
                                             q.boundary[1][ib0],
                                             q.boundary[1][ib1] );
            res[ic]+=fr;
        }
     }
      else if ( boundaryId==2 ) // bottom face
     {
         const size_t j=0;
         for( size_t i=0; i<ni; ++i )
        {
         // interior cell  index
            const CellIdx ic{i,j};

         // boundary cell index
            const CellIdx ib0{i,0};
            const CellIdx ib1{i,1};
   
         // face node indices
            const NodeIdx ip0{i+1,j};
            const NodeIdx ip1{i  ,j};
   
            const FluxRes fr = boundaryFlux( species,
                                             surface( mesh.nodes[ip0],
                                                      mesh.nodes[ip1] ),
                                             mesh.cells[ic],
                                             dq[ic][1],
                                             q.interior[ic],
                                             q.boundary[2][ib0],
                                             q.boundary[2][ib1] );
            res[ic]+=fr;
        }
     }
      else if ( boundaryId==3 ) // top face
     {
         const size_t j=nj-1;
         for( size_t i=0; i<ni; ++i )
        {
         // interior cell  index
            const CellIdx ic{i,j};

         // boundary cell index
            const CellIdx ib0{i,0};
            const CellIdx ib1{i,1};
   
         // face node indices
            const NodeIdx ip0{i  ,j+1};
            const NodeIdx ip1{i+1,j+1};
   
            const FluxRes fr = boundaryFlux( species,
                                             surface( mesh.nodes[ip0],
                                                      mesh.nodes[ip1] ),
                                             mesh.cells[ic],
                                             dq[ic][1],
                                             q.interior[ic],
                                             q.boundary[3][ib0],
                                             q.boundary[3][ib1] );
            res[ic]+=fr;
        }
     }
      else{ std::cout << boundaryId << "\n"; assert( false && "invalid boundary id for 2D, must be 0,1,2 or 3" ); }

      return;
  }

