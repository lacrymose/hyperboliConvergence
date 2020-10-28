
namespace geom
{
   template<floating_point Real>
   Surface1<Real> surface( const Point1<Real>& p )
  {
      return Surface1<Real>{.area   = Real(1.),
                            .centre = p,
                            .metric = Metric1<Real>{Direction1<Real>{Real(1.)}}};
  }

   template<floating_point Real>
   Surface2<Real> surface( const Point2<Real>& p0, const Point2<Real>& p1 )
  {
      using Mtrc = Metric2<Real>;
      using Dir  = Direction2<Real>;

      const Dir tangent0=p1-p0;
      const Dir tangent1=p0-p1;
      const Real l = length(tangent0);
      const Real l1= 1./l;

   // using orthog(tangent1) gives a right-handed system using fewer operations than using flip(orthog(tangent0))

      return {.area   = l,
              .centre = p0 + 0.5*tangent0,
              .metric = Mtrc{ orthog(l1*tangent1),
                                     l1*tangent0 }};
  }
}
