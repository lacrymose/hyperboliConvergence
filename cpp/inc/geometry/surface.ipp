
   Surface<1> surface( const Point<1>& p )
  {
      Surface<1> surf;
      surf.area   = 1;
      surf.centre = p;
      surf.metric = Metric<1>{Direction<1>{1.}};
      return surf;
  }

   Surface<2> surface( const Point<2>& p0, const Point<2>& p1 )
  {
      const Direction<2> tangent=p1-p0;

      Surface<2> surf;
      surf.area   = length( tangent );
      surf.centre = p0 + 0.5*tangent;
      surf.metric = Metric<2>{ cross(tangent),
                                     tangent };
      return surf;
  }
