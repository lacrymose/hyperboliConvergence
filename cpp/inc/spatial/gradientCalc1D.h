
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <cassert>

   template<ImplementedVarSet VarSetT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,1,Real,VarSetT>
   par::Array<vardelta_t<VarSetT>,1> gradientCalc( const par::Array<geom::Volume<1,Real>,1>& cells,
                                                   const par::Array<VarSetT,1>& q )
  {
      par::Array<vardelta_t<VarSetT>,1> dq(q.shape());
      gradientCalc( cells,q,dq );
      return dq;
  }

   template<ImplementedVarSet VarSetT, ImplementedVarDelta VarDelT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,1,Real,VarSetT,VarDelT>
   void gradientCalc( const par::Array<geom::Volume<1,Real>,1>& cells,
                      const par::Array<VarSetT,1>& q,
                            par::Array<VarDelT,1>& dq )
  {

      assert( cells.shape() ==  q.shape() );
      assert( cells.shape() == dq.shape() );
      const size_t nc = cells.shape(0);

      par::fill( dq, VarDelT{} );

   // zero derivative over boundaries
      for( size_t i=0; i<nc-1; i++ )
     {
         const VarDelT d = q[{i+1}] - q[{i}];

         dq[{i  }]+=d;
         dq[{i+1}]+=d;
     }

      return;
  }

/*
 * Periodic differences
 */
   template<ImplementedVarSet VarSetT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,1,Real,VarSetT>
   par::Array<vardelta_t<VarSetT>,1> gradientCalcP( const par::Array<geom::Volume<1,Real>,1>& cells,
                                                    const par::Array<VarSetT,1>& q )
  {
      par::Array<vardelta_t<VarSetT>,1> dq(q.shape());
      gradientCalcP( cells,q,dq );
      return dq;
  }

   template<ImplementedVarSet VarSetT, ImplementedVarDelta VarDelT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,1,Real,VarSetT,VarDelT>
   void gradientCalcP( const par::Array<geom::Volume<1,Real>,1>& cells,
                       const par::Array<VarSetT,1>& q,
                             par::Array<VarDelT,1>& dq )
  {

      assert( cells.shape() ==  q.shape() );
      assert( cells.shape() == dq.shape() );
      const size_t nc = cells.shape(0);

      par::fill( dq, VarDelT{} );

   // zero derivative over boundaries
      for( size_t i=0; i<nc-1; i++ )
     {
         const size_t il = i;
         const size_t ir = i+1;

         const VarDelT d = q[{ir}] - q[{il}];

         dq[{il}]+=d;
         dq[{ir}]+=d;
     }

   // periodic edge
     {
         const size_t il = nc-1;
         const size_t ir = 0;

         const VarDelT d = q[{ir}] - q[{il}];

         dq[{il}]+=d;
         dq[{ir}]+=d;
     }

      return;
  }

