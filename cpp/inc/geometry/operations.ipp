
# include <cmath>

namespace geom
{
   template<int nDim, floating_point Real>
   Real length2( const Direction<nDim,Real>& d )
  {
      Real l=0;
      for( const Real& x : d.elems ){ l+=x*x; }
      return l;
  }

   template<int nDim, floating_point Real>
   Real length( const Direction<nDim,Real>& d )
  {
      return std::sqrt(length2(d));
  }

   template<floating_point Real>
   Direction<2,Real> cross( const Direction<2,Real>& a )
  {
      return Direction<2,Real>{-a[1],a[0]};
  }

   template<floating_point Real>
   Direction<3,Real> cross( const Direction<3,Real>& a,
                            const Direction<3,Real>& b )
  {
      return Direction<3,Real>{ a[1]*b[2] - a[2]*b[1],
                                a[2]*b[0] - a[0]*b[2],
                                a[0]*b[1] - a[1]*b[0] };
  }
}
