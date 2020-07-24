
namespace geom
{
   template<floating_point Real>
   Surface<1,Real> surface( const Point<1,Real>& p )
  {
      return Surface<1,Real>{.area   = Real(1.),
                             .centre = p,
                             .metric = Metric<1,Real>{Direction<1,Real>{Real(1.)}}};
  }

   template<floating_point Real>
   Surface<2,Real> surface( const Point<2,Real>& p0, const Point<2,Real>& p1 )
  {
      const Direction<2,Real> tangent=p1-p0;
      const Real l = length(tangent);
      return Surface<2,Real>{.area   = l,
                             .centre = p0 + 0.5*tangent,
                             .metric = Metric<2,Real>{ -1.*norm(orthog(tangent)),
                                                        norm(       tangent ) }};
  }
}
