
namespace geom
{
/*
 * vertex centroid of a polygon is average of all vertices
 */
   template<typename Pt0, typename... Pts>
      requires    is_Point_v<Pt0>
              && ((std::is_same_v<Pt0,Pts>)&&...)
   Pt0 vertex_centroid( const Pt0& pt0,  const Pts&... pts )
  {
      constexpr int nDim = dim_of_v<Pt0>;
      using Real  = fptype_of_t<Pt0>;

      using Point = Point<    nDim,Real>;
      using Dir   = Direction<nDim,Real>;

      constexpr Real n1 = 1./( 1+sizeof...(Pts) );

      return Point{ n1*( Dir(pt0) + ((Dir(pts))+...) ) };
  }

/*
 * centre of mass of a quadrilateral
 *    make one line joining centroids of triangles from splitting along one   diagonal
 *    make one line joining centroids of triangles from splitting along other diagonal
 *    centre of mass is intersection of these two lines
 *
 *    2 ---- 3
 *    |      |
 *    |      |
 *    0 ---- 1
 *
 */
/*
   Is this an efficient way of calculating this? How about the FE basis element method from h4x?

   template<floating_point Real>
   Point<2,Real> mass_centroid( const Point<2,Real>& p0, const Point<2,Real>& p1,
                                const Point<2,Real>& p2, const Point<2,Real>& p3 )
  {
      using Point = Point<2,Real>;

      const Point p00 = vertex_centroid( p0,p1,p3 );
      const Point p01 = vertex_centroid( p0,p2,p3 );

      const Point p10 = vertex_centroid( p1,p0,p2 );
      const Point p11 = vertex_centroid( p1,p3,p2 );

      const std::array<Real,2> l0 = line_coefficients( p00, p01-p00 );
      const std::array<Real,2> l1 = line_coefficients( p10, p11-p10 );

      const Real x = ( l1[1] - l0[1] ) / ( l0[0] - l0[1] );
      const Real y = l0[0]*x + l0[1];

      return Point{x,y};
  }
*/
}
