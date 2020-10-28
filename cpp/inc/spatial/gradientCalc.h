
# pragma once

# include <solutionField/solutionField.h>
# include <conservationLaws/base/base.h>

# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <array>
# include <cassert>

   template<ImplementedVarSet VarSetT,
            int                  nDim,
            floating_point       Real>
      requires ConsistentTypes<law_of_v<VarSetT>,
                               nDim,
                               Real,
                               VarSetT>
   par::DualArray<vardelta_t<VarSetT>,nDim> gradientCalc( const Mesh<nDim,Real>&          mesh,
                                                          const SolutionField<VarSetT,nDim>& q )
  {
      par::DualArray<vardelta_t<VarSetT>,nDim> dq(q.shape());
      gradientCalc( mesh,q, dq );
      return dq;
  }

   template<ImplementedVarSet   VarSetT,
            ImplementedVarDelta VarDelT,
            int                    nDim,
            size_t                    N,
            floating_point         Real>
      requires ConsistentTypes<law_of_v<VarSetT>,
                               nDim,
                               Real,
                               VarSetT,
                               VarDelT>
               && N==nDim
   void gradientCalc( const Mesh<nDim,Real>&                          mesh,
                      const SolutionField<VarSetT,nDim>&                 q,
                            par::DualArray<std::array<VarDelT,N>,nDim>& dq )
  {
      par::fill( dq, std::array<VarDelT,N>{} );
      interiorGradient( mesh,q, dq );
      boundaryGradient( mesh,q, dq );
  }

/*
 * gradient calculation for interior of 1D domain
 */
   template<ImplementedVarSet   VarSetT,
            ImplementedVarDelta VarDelT,
            floating_point         Real>
      requires ConsistentTypes<law_of_v<VarSetT>,
                               1,
                               Real,
                               VarSetT,
                               VarDelT>
   void interiorGradient( const Mesh<1,Real>&                         mesh,
                          const SolutionField<VarSetT,1>&                q,
                                par::DualArray1<std::array<VarDelT,1>>& dq )
  {

      assert( mesh.cells.shape() == q.interior.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      const size_t nc = mesh.cells.shape(0);

      using CellIdx = typename SolutionField<VarSetT,1>::VarField::IdxType;

      for( size_t i=0; i<nc-1; ++i )
     {
         const CellIdx idx0{i};
         const CellIdx idx1{i+1};

         const VarDelT d = q.interior(idx1) - q.interior(idx0);

         dq(idx0)[0]+=d;
         dq(idx1)[0]+=d;
     }

      return;
  }

/*
 * gradient calculation for interior of 2D domain
 */
   template<ImplementedVarSet   VarSetT,
            ImplementedVarDelta VarDelT,
            floating_point         Real>
      requires ConsistentTypes<law_of_v<VarSetT>,
                               2,
                               Real,
                               VarSetT,
                               VarDelT>
   void interiorGradient( const Mesh<2,Real>&                         mesh,
                          const SolutionField<VarSetT,2>&                q,
                                par::DualArray2<std::array<VarDelT,2>>& dq )
  {
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

      par::fill( dq, std::array<VarDelT,2>{} );

      const size_t ni = mesh.cells.shape(0);
      const size_t nj = mesh.cells.shape(1);

      using CellIdx = typename SolutionField<VarSetT,2>::VarField::IdxType;

/*
      const auto dq_calc = []( const VarSetT& ql,
                               const VarSetT& qr ) -> VarDelT
     {
         return qr-ql;
     };

      const auto accumulator = []( const VarDelT& acc_old,
                                   const VarDelT& dq_new ) -> VarDelT
     {
         return acc_old+=dq_new;
     };

      par::accumulate_edge( parallel_schedule,
                            dq_calc,
                            accumulator,
                            accumulator,
                            dq,
                            q.interior );
*/

   // derivatives over i-normal faces
# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t j=0; j<nj; ++j )
     {
         for( size_t i=0; i<ni-1; ++i )
        {
            const CellIdx icl{i  ,j};
            const CellIdx icr{i+1,j};

            const VarDelT d = q.interior(icr) - q.interior(icl);
   
            dq(icl)[0]+=d;
            dq(icr)[0]+=d;
        }
     }

   // derivatives over i-normal faces
# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj-1; ++j )
        {
            const CellIdx icl{i,j  };
            const CellIdx icr{i,j+1};

            const VarDelT d = q.interior(icr) - q.interior(icl);
   
            dq(icl)[1]+=d;
            dq(icr)[1]+=d;
        }
     }

      return;
  }

/*
 * gradient calculation for boundary of 1D domain
 */
   template<ImplementedVarSet   VarSetT,
            ImplementedVarDelta VarDelT,
            floating_point         Real>
      requires ConsistentTypes<law_of_v<VarSetT>,
                               1,
                               Real,
                               VarSetT,
                               VarDelT>
   void boundaryGradient( const Mesh<1,Real>&                         mesh,
                          const SolutionField<VarSetT,1>&                q,
                                par::DualArray1<std::array<VarDelT,1>>& dq )
  {
      assert( mesh.cells.shape() == q.interior.shape() );
      assert( mesh.cells.shape() == dq.shape() );

      const size_t nc = mesh.cells.shape(0);

      using CellIdx = typename SolutionField<VarSetT,1>::VarField::IdxType;

   // non-zero gradients only over periodic boundary
      using BCType = BoundaryType<law_of_v<VarSetT>>;
      if( q.bcTypes[0] == BCType::Periodic )
     {
         assert( q.bcTypes[1] == BCType::Periodic );

         const CellIdx il{nc-1};
         const CellIdx ir{0};

         const VarDelT d = q.interior(ir) - q.interior(il);

         dq(il)[0]+=d;
         dq(ir)[0]+=d;
     }

      return;
  }

/*
 * gradient calculation for boundary of 2D domain
 */
   template<ImplementedVarSet   VarSetT,
            ImplementedVarDelta VarDelT,
            floating_point         Real>
      requires ConsistentTypes<law_of_v<VarSetT>,
                               2,
                               Real,
                               VarSetT,
                               VarDelT>
   void boundaryGradient( const Mesh<2,Real>&                         mesh,
                          const SolutionField<VarSetT,2>&                q,
                                par::DualArray2<std::array<VarDelT,2>>& dq )
  {
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

      const size_t ni = mesh.cells.shape(0);
      const size_t nj = mesh.cells.shape(1);

      using CellIdx = typename SolutionField<VarSetT,2>::VarField::IdxType;

   // non-zero gradients only over periodic boundary
      using BCType = BoundaryType<law_of_v<VarSetT>>;

      if( q.bcTypes[0] == BCType::Periodic ) // periodic boundary along +/- i-direction
     {
      // check consistency
         assert( q.bcTypes[1] == BCType::Periodic );

         for( size_t j=0; j<nj-1; ++j )
        {
            const size_t il = ni-1;
            const size_t ir = 0;
   
         // cell left/right indices
            const CellIdx icl{il,j};
            const CellIdx icr{ir,j};
   
            const VarDelT d = q.interior(icr) - q.interior(icl);
   
            dq(icl)[0]+=d;
            dq(icr)[0]+=d;
        }
     }

      if( q.bcTypes[2] == BCType::Periodic ) // periodic boundary along +/- j-direction
     {
      // check consistency
         assert( q.bcTypes[3] == BCType::Periodic );

         for( size_t i=0; i<ni-1; ++i )
        {
            const size_t jl = nj-1;
            const size_t jr = 0;
   
         // cell left/right indices
            const CellIdx icl{i,jl};
            const CellIdx icr{i,jr};
   
            const VarDelT d = q.interior(icr) - q.interior(icl);
   
            dq(icl)[1]+=d;
            dq(icr)[1]+=d;
        }
     }

      return;
  }
