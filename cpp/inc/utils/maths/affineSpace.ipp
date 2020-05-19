
// --------------- Point in-place arithmetic ---------------

   // p+=d
   template<int NDIM, typename Point, typename Delta, floating_point Real>
   Point& AffinePointBase<NDIM,Point,Delta,Real>::operator+=( const Delta& d )
  {
      for( int i=0; i<NDIM; i++ ){ v[i]+=d[i]; }
      return static_cast<Point&>(*this);
  }

   // p-=d
   template<int NDIM, typename Point, typename Delta, floating_point Real>
   Point& AffinePointBase<NDIM,Point,Delta,Real>::operator-=( const Delta& d )
  {
      for( int i=0; i<NDIM; i++ ){ v[i]-=d[i]; }
      return static_cast<Point&>(*this);
  }


// --------------- Delta in-place arithmetic ---------------

   // d+=d
   template<int NDIM, typename Point, typename Delta, floating_point Real>
   Delta& AffineDeltaBase<NDIM,Point,Delta,Real>::operator+=( const Delta& d )
  {
      for( int i=0; i<NDIM; i++ ){ v[i]+=d[i]; }
      return static_cast<Delta&>(*this);
  }

   // d-=d
   template<int NDIM, typename Point, typename Delta, floating_point Real>
   Delta& AffineDeltaBase<NDIM,Point,Delta,Real>::operator-=( const Delta& d )
  {
      for( int i=0; i<NDIM; i++ ){ v[i]-=d[i]; }
      return static_cast<Delta&>(*this);
  }

   // d*=a
   template<int NDIM, typename Point, typename Delta, floating_point Real>
   Delta& AffineDeltaBase<NDIM,Point,Delta,Real>::operator*=( const Real a )
  {
      for( Real& w : v ){ w*=a; }
      return static_cast<Delta&>(*this);
  }

   // d/=a
   template<int NDIM, typename Point, typename Delta, floating_point Real>
   Delta& AffineDeltaBase<NDIM,Point,Delta,Real>::operator/=( const Real a )
  {
      const Real a1=1./a;
      for( Real& w : v ){ w*=a1; }
      return static_cast<Delta&>(*this);
  }


// --------------- Point arithmetic ---------------

   // p = p+d
   template<typename Point, typename Delta>
      requires is_affine_pair<Point,Delta>::value
   Point operator+( const Point& p, const Delta& d )
  {
      Point result(p);
      result+=d;
      return result;
  }

   // p = p-d
   template<typename Point, typename Delta>
      requires is_affine_pair<Point,Delta>::value
   Point operator-( const Point& p, const Delta& d )
  {
      Point result(p);
      result-=d;
      return result;
  }


// --------------- Delta arithmetic ---------------

   // d = d+d
   template<typename Delta>
      requires has_affinedelta_base<Delta>::value
   Delta operator+( const Delta& lhs, const Delta& rhs )
  {
      Delta result(lhs);
      result+=rhs;
      return result;
  }

   // d = d-d
   template<typename Delta>
      requires has_affinedelta_base<Delta>::value
   Delta operator-( const Delta& lhs, const Delta& rhs )
  {
      Delta result(lhs);
      result-=rhs;
      return result;
  }

   // d = p-p
   template<typename Point>
      requires has_affinepoint_base<Point>::value
   typename Point::delta_type operator-( const Point& lhs, const Point& rhs )
  {
      using Delta = typename Point::delta_type;
      Delta result;
      for( int i=0; i<Point::N; i++ ){ result[i]=lhs[i]-rhs[i]; }
      return result;
  }

   // d = a*d
   template<typename Delta>
      requires has_affinedelta_base<Delta>::value
   Delta operator*( const typename Delta::value_type a, const Delta& d )
  {
      Delta result(d);
      result*=a;
      return result;
  }

   // d = d*a
   template<typename Delta>
      requires has_affinedelta_base<Delta>::value
   Delta operator*( const Delta& d, const typename Delta::value_type a )
  {
      Delta result(d);
      result*=a;
      return result;
  }

   // d = d/a
   template<typename Delta, floating_point Real>
      requires has_affinedelta_base<Delta>::value
   Delta operator/( const Delta& d, const typename Delta::value_type a )
  {
      const typename Delta::value_type a1=1./a;
      Delta result(d);
      result*=a1;
      return result;
  }


// --------------- Printing to stream ---------------

/*
 * send each element of Point to stream, seperated by a single space
 */
   template<typename Point>
      requires has_affinepoint_base<Point>::value
   std::ostream& operator<<( std::ostream& os, const Point& p )
  {
      constexpr int N=Point::N;

      for( int i=0; i<N-1; i++ )
     {
         os << p[i] << " ";
     }
      os << p[N-1];
      return os;
  }

/*
 * send each element of Direction to stream, seperated by a single space
 */
   template<typename Delta>
      requires has_affinedelta_base<Delta>::value
   std::ostream& operator<<( std::ostream& os, const Delta& d )
  {
      constexpr int N=Delta::N;

      for( int i=0; i<N-1; i++ )
     {
         os << d[i] << " ";
     }
      os << d[N-1];
      return os;
  }


