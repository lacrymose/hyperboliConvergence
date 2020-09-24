
namespace geom
{
   template<floating_point Real>
   Volume<1,Real> volume( const Point<1,Real>& p0, const Point<1,Real>& p1 )
  {
      const Direction<1,Real> dir=p1-p0;
      return Volume<1,Real>{.volume=length(dir),
                            .centre=p0+0.5*dir};
  }

   // corners of standard quad
   //
   //   2 ---- 3
   //   |      |
   //   |      |
   //   0 ---- 1
   //
   template<floating_point Real>
   Volume<2,Real> volume( const Point<2,Real>& p0, const Point<2,Real>& p1,
                          const Point<2,Real>& p2, const Point<2,Real>& p3 )
  {
      using Dir = Direction<2,Real>;

   // distance vectors along each side
      const Dir d01 = p1 - p0;
      const Dir d02 = p2 - p0;

      const Dir d31 = p1 - p3;
      const Dir d32 = p2 - p3;
 
      return Volume<2,Real>{.volume = Real(0.25*cross(d01,d02)*cross(d32,d31)),
                            .centre = vertex_centroid( p0,p1,p2,p3 ) };
// this should be mass centroid but don't have efficient implementation yet. vertex centroid ok for near-orthogonal layouts
                         /* .centre = mass_centroid( p0,p1,p2,p3 ) }; */
  }
}
