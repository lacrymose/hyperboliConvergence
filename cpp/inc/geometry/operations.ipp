
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
   Real cross( const Direction2<Real>& a,
               const Direction2<Real>& b )
  {
      return a[0]*b[1] - a[1]*b[0];
  }

   template<floating_point Real>
   Direction3<Real> cross( const Direction3<Real>& a,
                           const Direction3<Real>& b )
  {
      return Direction3<Real>{ a[1]*b[2] - a[2]*b[1],
                               a[2]*b[0] - a[0]*b[2],
                               a[0]*b[1] - a[1]*b[0] };
  }

   template<floating_point Real>
   Direction2<Real> orthog( const Direction2<Real>& d )
  {
      return {-d[1],d[0]};
  }

   template<int nDim, floating_point Real>
   Direction<nDim,Real> norm( const Direction<nDim,Real>& d )
  {
      return d/length(d);
  }

   template<floating_point Real>
   std::array<Real,2> line_coefficients( const Point2<Real>& p, const Direction2<Real>& d )
  {
      const Real a = d[1]/d[0];
      return {a, p[1] - a*p[0]};
  }

   template<int nDim, floating_point Real>
   Direction<nDim,Real> flip ( const Direction<nDim,Real>& d0 )
  {
      return -1.*d0;
  }

   template<int nDim, floating_point Real>
   Surface<nDim,Real> flip ( const Surface<nDim,Real>& s0 )
  {
      Surface<nDim,Real> s1(s0);
      if constexpr( nDim==1 )
     {
         s1.metric[0] = flip( s0.metric[0] );
     }
      else if constexpr( nDim==2 )
     {
         s1.metric[0] = flip( s0.metric[0] );
         s1.metric[1] = flip( s0.metric[1] );
     }
      else if constexpr( nDim==3 )
     {
         s1.metric[0] = flip( s0.metric[0] );
         s1.metric[1] = flip( s0.metric[1] );
     }
      else
     {
         constexpr bool validDim = nDim==1 or nDim==2 or nDim==3;
         static_assert( validDim, "flipping Surface<nDim,Real> with nDim>3 not defined" );
     }
      return s1;
  }
}
