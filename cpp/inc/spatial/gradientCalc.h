
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <mdarray/mdarray.h>

# include <cassert>

   template<ImplementedVarSet VarSetT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,1,Real,VarSetT>
   MDArray<vardelta_t<VarSetT>,1> gradientCalc( const MDArray<geom::Volume<1,Real>,1>& cells,
                                                const MDArray<VarSetT,1>& q )
  {
      MDArray<vardelta_t<VarSetT>,1> dq(q.dims);
      gradientCalc( cells,q,dq );
      return dq;
  }

   template<ImplementedVarSet VarSetT, ImplementedVarDelta VarDelT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,1,Real,VarSetT,VarDelT>
   void gradientCalc( const MDArray<geom::Volume<1,Real>,1>& cells,
                      const MDArray<VarSetT,1>& q,
                            MDArray<VarDelT,1>& dq )
  {

      assert( cells.dims ==  q.dims );
      assert( cells.dims == dq.dims );
      const size_t nc = cells.dims[0];

      for( VarDelT& d : dq.elems ){ d = VarDelT{}; }

   // zero derivative over boundaries
      for( size_t i=0; i<nc-1; i++ )
     {
         const VarDelT d = q[{i+1}] - q[{i}];

         dq[{i  }]+=d;
         dq[{i+1}]+=d;
     }

      return;
  }
