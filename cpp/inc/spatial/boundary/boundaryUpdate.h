
# pragma once

# include <spatial/boundary/boundaryCondition.h>

# include <solutionField/solutionField.h>
# include <conservationLaws/base/base.h>

# include <mesh/mesh.h>
# include <geometry/geometry.h>

# include <utils/utils.h>
# include <utils/concepts.h>

# include <tuple>
# include <cassert>


   template<LawType                  Law,
            int                     nDim,
            ImplementedVarSet    SolVarT,
            typename...    BoundaryConds,
            floating_point          Real>
      requires ConsistentTypes<Law,
                               nDim,
                               Real,
                               SolVarT>
   void boundaryUpdate( const Mesh<nDim,Real>&             mesh,
                        const std::tuple<BoundaryConds...>  bcs,
                        const Species<Law,Real>&        species,
                              SolutionField<SolVarT,nDim>&    q )
  {
   // check mesh sizes
      assert( mesh.cells.shape() == q.interior.shape() );

      const auto call_bc_update = []( auto&&... args )
     {
         boundaryUpdate( std::forward<decltype(args)>(args)... );
     };

   // for each boundary, update the values from the matching boundary condition in the bc tuple
      for( unsigned int boundaryId=0; boundaryId<q.nBoundaries; ++boundaryId )
     {
         selectBoundaryCondition( q.bcTypes[boundaryId], bcs, call_bc_update,
                                  boundaryId, mesh, species, q );
     }
      return;
  }

/*
 * boundary values for periodic boundaries are unused
 */
   template<LawType               Law,
            int                  nDim,
            ImplementedVarSet SolVarT,
            floating_point       Real>
      requires ConsistentTypes<Law,
                               nDim,
                               Real,
                               SolVarT>
   void boundaryUpdate( const BoundaryCondition<Law,BoundaryType<Law>::Periodic>&,
                        const size_t              boundaryId,
                        const Mesh<nDim,Real>&          mesh,
                        const Species<Law,Real>&     species,
                              SolutionField<SolVarT,nDim>& q )
  {
      return;
  }

/*
 * boundary values for 1D non-periodic boundaries are updated using function given by boundary condition
 */
   template<LawType                   Law,
            BoundaryType<Law>      BCType,
            ImplementedVarSet     SolVarT,
            typename           UpdateFunc,
            typename...      OtherBCFuncs,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 1,
                                 Real,
                                 SolVarT>
   void boundaryUpdate( const BoundaryCondition<Law,
                                                BCType,
                                                UpdateFunc,
                                                OtherBCFuncs...>& bc,
                        const size_t                      boundaryId,
                        const Mesh<1,Real>&                     mesh,
                        const Species<Law,Real>&             species,
                              SolutionField<SolVarT,1>&            q )
  {
   // check array sizes
      assert( mesh.cells.shape() == q.interior.shape() );

   // check boundary condition selection was correct
      assert( q.bcTypes[boundaryId] == BCType );

      const auto& updateBC = std::get<0>(bc.funcs);

      using CellIdx = typename SolutionField<SolVarT,1>::VarField::IdxType;
      using NodeIdx = typename Mesh<1,Real>::NodeArray::IdxType;

      const size_t bID = boundaryId;
      if( boundaryId==0 ) // left boundary
     {
         const NodeIdx ip{0};    // node on boundary
         const CellIdx ic{0};    // interior cell
         const CellIdx ibv{0};   // boundary value
         const CellIdx ibr{1};   // boundary reference value

         q.boundary[bID](ibv) = updateBC( species,
                                          surface( mesh.nodes(ip) ), mesh.cells(ic),
                                          q.interior(ic), q.boundary[bID](ibr) );
     }
      else if ( boundaryId==1 ) // right boundary
     {
         const size_t nc = q.interior.shape(0);

         const NodeIdx ip{nc};     // node on boundary
         const CellIdx ic{nc-1};   // interior cell
         const CellIdx ibv{0};     // boundary value
         const CellIdx ibr{1};     // boundary reference value

      // boundary surface must face into domain
         q.boundary[bID](ibv) = updateBC( species,
                                          flip( surface( mesh.nodes(ip) ) ), mesh.cells(ic),
                                          q.interior(ic), q.boundary[bID](ibr) );
     }
      else{ assert( false && "invalid boundary id for 1D boundary update, must be 0 or 1" ); }

      return;
  }

/*
 * boundary values for 2D non-periodic boundaries are updated using function given by boundary condition
 */
   template<LawType                   Law,
            BoundaryType<Law>      BCType,
            ImplementedVarSet     SolVarT,
            typename           UpdateFunc,
            typename...      OtherBCFuncs,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 2,
                                 Real,
                                 SolVarT>
   void boundaryUpdate( const BoundaryCondition<Law,
                                                BCType,
                                                UpdateFunc,
                                                OtherBCFuncs...>& bc,
                        const size_t                      boundaryId,
                        const Mesh<2,Real>&                     mesh,
                        const Species<Law,Real>&             species,
                              SolutionField<SolVarT,2>&            q )
  {
   // check array sizes
      assert( mesh.cells.shape() == q.interior.shape() );

   // check boundary condition selection was correct
      assert( q.bcTypes[boundaryId] == BCType );

      const auto& updateBC = std::get<0>(bc.funcs);

      const size_t ni = mesh.cells.shape(0);
      const size_t nj = mesh.cells.shape(1);

      using CellIdx = typename SolutionField<SolVarT,2>::VarField::IdxType;
      using NodeIdx = typename Mesh<2,Real>::NodeArray::IdxType;

      const size_t bID = boundaryId;
      if( boundaryId==0 ) // left boundary
     {
         const size_t i=0;
         for( size_t j=0; j<nj; ++j )
        {
         // interior cell  index
            const CellIdx ic{i,j};

         // boundary cell index
            const CellIdx ibv{j,0};
            const CellIdx ibr{j,1};
   
         // face node indices
            const NodeIdx ip0{i,j  };
            const NodeIdx ip1{i,j+1};
   
            q.boundary[bID](ibv) = updateBC( species,
                                             surface( mesh.nodes(ip0),
                                                      mesh.nodes(ip1) ),
                                             mesh.cells(ic),
                                             q.interior(ic), q.boundary[bID](ibr) );
        }
     }
      else if ( boundaryId==1 ) // right boundary
     {
         const size_t i=ni-1;
         for( size_t j=0; j<nj; ++j )
        {
         // interior cell  index
            const CellIdx ic{i,j};

         // boundary cell index
            const CellIdx ibv{j,0};
            const CellIdx ibr{j,1};
   
         // face node indices
            const NodeIdx ip0{i+1,j+1};
            const NodeIdx ip1{i+1,j  };
   
            q.boundary[bID](ibv) = updateBC( species,
                                             surface( mesh.nodes(ip0),
                                                      mesh.nodes(ip1) ),
                                             mesh.cells(ic),
                                             q.interior(ic), q.boundary[bID](ibr) );
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
            const CellIdx ibv{i,0};
            const CellIdx ibr{i,1};
   
         // face node indices
            const NodeIdx ip0{i+1,j};
            const NodeIdx ip1{i  ,j};
   
            q.boundary[bID](ibv) = updateBC( species,
                                             surface( mesh.nodes(ip0),
                                                      mesh.nodes(ip1) ),
                                             mesh.cells(ic),
                                             q.interior(ic), q.boundary[bID](ibr) );
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
            const CellIdx ibv{i,0};
            const CellIdx ibr{i,1};
   
         // face node indices
            const NodeIdx ip0{i  ,j+1};
            const NodeIdx ip1{i+1,j+1};
   
            q.boundary[bID](ibv) = updateBC( species,
                                             surface( mesh.nodes(ip0),
                                                      mesh.nodes(ip1) ),
                                             mesh.cells(ic),
                                             q.interior(ic), q.boundary[bID](ibr) );
        }
     }
      else{ assert( false && "invalid boundary id for 2D boundary update, must be 0,1,2 or 3" ); }

      return;
  }

