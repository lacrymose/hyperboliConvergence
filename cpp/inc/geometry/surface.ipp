
namespace Geometry
{
   template<floating_point Real>
   Surface<1,Real> surface( const Point<1,Real>& p )
  {
      Surface<1,Real> surf;
      surf.area   = 1;
      surf.centre = p;
      surf.metric = Metric<1,Real>{Direction<1,Real>{1.}};
      return surf;
  }

   template<floating_point Real>
   Surface<2,Real> surface( const Point<2,Real>& p0, const Point<2,Real>& p1 )
  {
      const Direction<2,Real> tangent=p1-p0;

      Surface<2,Real> surf;
      surf.area   = length( tangent );
      surf.centre = p0 + 0.5*tangent;
      surf.metric = Metric<2,Real>{ cross(tangent),
                                          tangent };
      return surf;
  }
}
