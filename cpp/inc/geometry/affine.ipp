
/*
 * Point<nDim> and Direction<nDim> form an affine space.
 *    8 arithmetic operations need defining for the affine space
 *    directions are made by subtracting points
 *    direction can be added/subtracted from other directions or points
 *    directions can be multiplied/divided by a scalar
 */

// d = d+d
   template<int nDim>
   Direction<nDim> operator+( const Direction<nDim>& d0, const Direction<nDim>& d1 )
  {
      Direction<nDim> d2(d0);
      d2+=d1;
      return d2;
  }

// d = d-d
   template<int nDim>
   Direction<nDim> operator-( const Direction<nDim>& d0, const Direction<nDim>& d1 )
  {
      Direction<nDim> d2(d0);
      d2-=d1;
      return d2;
  }

// d = p-p
   template<int nDim>
   Direction<nDim> operator-( const Point<nDim>& p0, const Point<nDim>& p1 )
  {
      using Dir = Direction<nDim>;
      return Dir(p0) - Dir(p1);
  }

// p = p+d
   template<int nDim>
   Point<nDim> operator+( const Point<nDim>& p0, const Direction<nDim>& d0 )
  {
      Point<nDim> p1(p0);
      p1+=d0;
      return p1;
  }

// p = p-d
   template<int nDim>
   Point<nDim> operator-( const Point<nDim>& p0, const Direction<nDim>& d0 )
  {
      Point<nDim> p1(p0);
      p1-=d0;
      return p1;
  }

// d = d*a
   template<int nDim>
   Direction<nDim> operator*( const Direction<nDim>& d0, const Types::Real a )
  {
      Direction<nDim> d1(d0);
      d1*=a;
      return d1;
  }

// d = a*d
   template<int nDim>
   Direction<nDim> operator*( const Types::Real a, const Direction<nDim>& d0 )
  {
      Direction<nDim> d1(d0);
      d1*=a;
      return d1;
  }

// d = d/a
   template<int nDim>
   Direction<nDim> operator/( const Direction<nDim>& d0, const Types::Real a )
  {
      Direction<nDim> d1(d0);
      d1/=a;
      return d1;
  }
