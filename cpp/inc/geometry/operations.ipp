
# include <cmath>

namespace geom
{
   template<int nDim, floating_point Real>
   Real dot( const Direction<nDim,Real>& l,
             const Direction<nDim,Real>& r )
  {
      Real d=0;
      for( int i=0; i<nDim; i++ ){ d+=l[i]*r[i]; }
      return d;
  }

   template<int nDim, floating_point Real>
   Real length2( const Direction<nDim,Real>& d )
  {
      return dot( d,d );
  }

   template<int nDim, floating_point Real>
   Real length( const Direction<nDim,Real>& d )
  {
      return std::sqrt(length2(d));
  }

   template<floating_point Real>
   Real cross( const Direction<2,Real>& a,
               const Direction<2,Real>& b )
  {
      return a[0]*b[1] - a[1]*b[0];
  }

   template<floating_point Real>
   Direction<3,Real> cross( const Direction<3,Real>& a,
                            const Direction<3,Real>& b )
  {
      return Direction<3,Real>{ a[1]*b[2] - a[2]*b[1],
                                a[2]*b[0] - a[0]*b[2],
                                a[0]*b[1] - a[1]*b[0] };
  }

   template<floating_point Real>
   Direction<2,Real> orthog( const Direction<2,Real>& d )
  {
      return Direction<2,Real>{-d[1],d[0]};
  }

   template<int nDim, floating_point Real>
   Direction<nDim,Real> norm( const Direction<nDim,Real>& d )
  {
      return d/length(d);
  }

   template<floating_point Real>
   std::array<Real,2> line_coefficients( const Point<2,Real>& p, const Direction<2,Real>& d )
  {
      const Real a = d[1]/d[0];
      return {a, p[1] - a*p[0]};
  }
}
