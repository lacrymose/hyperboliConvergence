
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <cassert>

   template<ImplementedVarSet VarSetT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,2,Real,VarSetT>
   par::Array<vardelta_t<VarSetT>,2> gradientCalc( const par::Array<geom::Volume<2,Real>,2>& cells,
                                                   const par::Array<VarSetT,2>& q )
  {
      par::Array<std::array<vardelta_t<VarSetT>,2>,2> dq(q.shape());
      gradientCalc( cells,q,dq );
      return dq;
  }

   template<ImplementedVarSet VarSetT, ImplementedVarDelta VarDelT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,2,Real,VarSetT,VarDelT>
   void gradientCalc( const par::Array<geom::Volume<2,Real>,2>& cells,
                      const par::Array<VarSetT,2>&                  q,
                            par::Array<std::array<VarDelT,2>,2>&   dq )
  {

      assert( cells.shape() ==  q.shape() );
      assert( cells.shape() == dq.shape() );

      par::fill( dq, std::array<VarDelT,2>{} );

      const size_t ni = cells.shape(0);
      const size_t nj = cells.shape(1);

   // zero derivative over boundaries

   // derivatives over i-normal faces
      for( size_t i=0; i<ni-1; i++ )
     {
         for( size_t j=0; j<nj; j++ )
        {
            const par::Idx<2> icl{i  ,j};
            const par::Idx<2> icr{i+1,j};

            const VarDelT d = q[icr] - q[icl];
   
            dq[icl][0]+=d;
            dq[icr][0]+=d;
        }
     }

   // derivatives over i-normal faces
      for( size_t i=0; i<ni; i++ )
     {
         for( size_t j=0; j<nj-1; j++ )
        {
            const par::Idx<2> icl{i,j  };
            const par::Idx<2> icr{i,j+1};

            const VarDelT d = q[icr] - q[icl];
   
            dq[icl][1]+=d;
            dq[icr][1]+=d;
        }
     }

      return;
  }

/*
 * Periodic differences
 */
   template<ImplementedVarSet VarSetT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,2,Real,VarSetT>
   par::Array<std::array<vardelta_t<VarSetT>,2>,2>
      gradientCalcP( const par::Array<geom::Volume<2,Real>,2>& cells,
                     const par::Array<VarSetT,2>& q )
  {
      par::Array<std::array<vardelta_t<VarSetT>,2>,2> dq(q.shape());
      gradientCalcP( cells,q,dq );
      return dq;
  }

   template<ImplementedVarSet VarSetT, ImplementedVarDelta VarDelT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,2,Real,VarSetT,VarDelT>
   void gradientCalcP( const par::Array<geom::Volume<2,Real>,2>& cells,
                       const par::Array<VarSetT,2>&                  q,
                             par::Array<std::array<VarDelT,2>,2>&   dq )
  {

      assert( cells.shape() ==  q.shape() );
      assert( cells.shape() == dq.shape() );

   // calculate gradient over internal faces
      gradientCalc( cells, q, dq );

      const size_t ni = cells.shape(0);
      const size_t nj = cells.shape(1);


   // periodic boundary along +/- i-direction
      for( size_t j=0; j<nj-1; j++ )
     {
         const size_t il = ni-1;
         const size_t ir = 0;

      // cell left/right indices
         const par::Idx<2> icl{il,j};
         const par::Idx<2> icr{ir,j};

         const VarDelT d = q[icr] - q[icl];

         dq[icl][0]+=d;
         dq[icr][0]+=d;
     }

   // periodic boundary along +/- j-direction
      for( size_t i=0; i<ni-1; i++ )
     {
         const size_t jl = nj-1;
         const size_t jr = 0;

      // cell left/right indices
         const par::Idx<2> icl{i,jl};
         const par::Idx<2> icr{i,jr};

         const VarDelT d = q[icr] - q[icl];

         dq[icl][1]+=d;
         dq[icr][1]+=d;
     }

      return;
  }

