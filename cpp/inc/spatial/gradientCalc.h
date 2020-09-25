
# pragma once

# include <solutionField/solutionField.h>
# include <conservationLaws/base/base.h>

# include <mesh/mesh.h>
# include <geometry/geometry.h>

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
   par::Array<vardelta_t<VarSetT>,nDim> gradientCalc( const Mesh<nDim,Real>&          mesh,
                                                      const SolutionField<VarSetT,nDim>& q )
  {
      par::Array<vardelta_t<VarSetT>,nDim> dq(q.shape());
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
   void gradientCalc( const Mesh<nDim,Real>&           mesh,
                      const SolutionField<VarSetT,nDim>&  q,
                            par::Array<std::array<VarDelT,N>,nDim>& dq )
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
   void interiorGradient( const Mesh<1,Real>&           mesh,
                          const SolutionField<VarSetT,1>&  q,
                                par::Array<std::array<VarDelT,1>,1>& dq )
  {

      assert( mesh.cells.shape() == q.interior.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      const size_t nc = mesh.cells.shape(0);

      for( size_t i=0; i<nc-1; i++ )
     {
         const VarDelT d = q.interior[{i+1}] - q.interior[{i}];
         dq[{i  }][0]+=d;
         dq[{i+1}][0]+=d;
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
   void interiorGradient( const Mesh<2,Real>&                      mesh,
                          const SolutionField<VarSetT,2>&             q,
                                par::Array<std::array<VarDelT,2>,2>& dq )
  {
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

      par::fill( dq, std::array<VarDelT,2>{} );

      const size_t ni = mesh.cells.shape(0);
      const size_t nj = mesh.cells.shape(1);

   // derivatives over i-normal faces
# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t j=0; j<nj; j++ )
     {
         for( size_t i=0; i<ni-1; i++ )
        {
            const par::Idx<2> icl{i  ,j};
            const par::Idx<2> icr{i+1,j};

            const VarDelT d = q.interior[icr] - q.interior[icl];
   
            dq[icl][0]+=d;
            dq[icr][0]+=d;
        }
     }

   // derivatives over i-normal faces
# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni; i++ )
     {
         for( size_t j=0; j<nj-1; j++ )
        {
            const par::Idx<2> icl{i,j  };
            const par::Idx<2> icr{i,j+1};

            const VarDelT d = q.interior[icr] - q.interior[icl];
   
            dq[icl][1]+=d;
            dq[icr][1]+=d;
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
   void boundaryGradient( const Mesh<1,Real>&           mesh,
                          const SolutionField<VarSetT,1>&  q,
                                par::Array<std::array<VarDelT,1>,1>& dq )
  {
      assert( mesh.cells.shape() == q.interior.shape() );
      assert( mesh.cells.shape() == dq.shape() );

      const size_t nc = mesh.cells.shape(0);

   // non-zero gradients only over periodic boundary
      using BCType = BoundaryType<law_of_v<VarSetT>>;
      if( q.bcTypes[0] == BCType::Periodic )
     {
         assert( q.bcTypes[1] == BCType::Periodic );

         const size_t il = nc-1;
         const size_t ir = 0;

         const VarDelT d = q.interior[{ir}] - q.interior[{il}];

         dq[{il}][0]+=d;
         dq[{ir}][0]+=d;
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
   void boundaryGradient( const Mesh<2,Real>&                      mesh,
                          const SolutionField<VarSetT,2>&             q,
                                par::Array<std::array<VarDelT,2>,2>& dq )
  {
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

      const size_t ni = mesh.cells.shape(0);
      const size_t nj = mesh.cells.shape(1);

   // non-zero gradients only over periodic boundary
      using BCType = BoundaryType<law_of_v<VarSetT>>;

      if( q.bcTypes[0] == BCType::Periodic ) // periodic boundary along +/- i-direction
     {
      // check consistency
         assert( q.bcTypes[1] == BCType::Periodic );

         for( size_t j=0; j<nj-1; j++ )
        {
            const size_t il = ni-1;
            const size_t ir = 0;
   
         // cell left/right indices
            const par::Idx<2> icl{il,j};
            const par::Idx<2> icr{ir,j};
   
            const VarDelT d = q.interior[icr] - q.interior[icl];
   
            dq[icl][0]+=d;
            dq[icr][0]+=d;
        }
     }

      if( q.bcTypes[2] == BCType::Periodic ) // periodic boundary along +/- j-direction
     {
      // check consistency
         assert( q.bcTypes[3] == BCType::Periodic );

         for( size_t i=0; i<ni-1; i++ )
        {
            const size_t jl = nj-1;
            const size_t jr = 0;
   
         // cell left/right indices
            const par::Idx<2> icl{i,jl};
            const par::Idx<2> icr{i,jr};
   
            const VarDelT d = q.interior[icr] - q.interior[icl];
   
            dq[icl][1]+=d;
            dq[icr][1]+=d;
        }
     }

      return;
  }
