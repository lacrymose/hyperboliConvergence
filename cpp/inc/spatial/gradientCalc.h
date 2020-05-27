
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <vector>

# include <cassert>

   template<ImplementedVarSet VarSetT, ImplementedVarDelta VarDelT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,1,Real,VarSetT>
   std::vector<vardelta_t<VarSetT>> gradientCalc( const std::vector<geom::Volume<1,Real>>& cells,
                                                  const std::vector<VarSetT>& q )
  {
      std::vector<vardelta_t<VarSetT>> dq(q.size());
      gradientCalc( cells,q,dq );
      return dq;
  }

   template<ImplementedVarSet VarSetT, ImplementedVarDelta VarDelT, floating_point Real>
      requires ConsistentTypes<law_of_v<VarSetT>,1,Real,VarSetT,VarDelT>
   void gradientCalc( const std::vector<geom::Volume<1,Real>>& cells,
                      const std::vector<VarSetT>& q,
                            std::vector<VarDelT>& dq )
  {
      const size_t nc = cells.size();

      assert( nc ==  q.size() );
      assert( nc == dq.size() );

      for( VarDelT& d : dq ){ d = VarDelT{}; }

   // zero derivative over boundaries
      for( size_t i=0; i<nc-1; i++ )
     {
         const VarDelT d = q[i+1] - q[i];

         dq[i  ]+=d;
         dq[i+1]+=d;
     }

      return;
  }
