
// --------------- Vector in-place arithmetic ---------------

   // v+=v
   template<int N, typename Derived, floating_point Real>
   Derived& VectorSpaceBase<N,Derived,Real>::operator+=( const Derived& other )
  {
      for( int i=0; i<N; i++ ){ v[i]+=other[i]; }
      return static_cast<Derived&>(*this);
  };

   // v-=v
   template<int N, typename Derived, floating_point Real>
   Derived& VectorSpaceBase<N,Derived,Real>::operator-=( const Derived& other )
  {
      for( int i=0; i<N; i++ ){ v[i]-=other[i]; }
      return static_cast<Derived&>(*this);
  };

   // v*=a
   template<int N, typename Derived, floating_point Real>
   Derived& VectorSpaceBase<N,Derived,Real>::operator*=( const Real a )
  {
      for( Real& w : v ){ w*=a; }
      return static_cast<Derived&>(*this);
  };

   // v/=a
   template<int N, typename Derived, floating_point Real>
   Derived& VectorSpaceBase<N,Derived,Real>::operator/=( const Real a )
  {
      const Real a1=1./a;
      for( Real& w : v ){ w*=a1; }
      return static_cast<Derived&>(*this);
  };


// --------------- Vector arithmetic ---------------

   // v = v+v
   template<typename Vec>
      requires has_vectorspace_base<Vec>::value
   Vec operator+( const Vec& lhs, const Vec& rhs )
  {
      Vec result(lhs);
      result+=rhs;
      return result;
  }

   // v = v-v
   template<typename Vec>
      requires has_vectorspace_base<Vec>::value
   Vec operator-( const Vec& lhs, const Vec& rhs )
  {
      Vec result(lhs);
      result-=rhs;
      return result;
  }

   // v = a*v
   template<typename Vec>
      requires has_vectorspace_base<Vec>::value
   Vec operator*( const typename Vec::value_type a, const Vec& rhs )
  {
      Vec result(rhs);
      result*=a;
      return result;
  }

   // v = v*a
   template<typename Vec>
      requires has_vectorspace_base<Vec>::value
   Vec operator*( const Vec& lhs, const typename Vec::value_type a )
  {
      Vec result(lhs);
      result*=a;
      return result;
  }

   // v = v/a
   template<typename Vec>
      requires has_vectorspace_base<Vec>::value
   Vec operator/( const Vec& lhs, const typename Vec::value_type a )
  {
      const typename Vec::value_type a1 = 1./a;
      Vec result(lhs);
      result*=a1;
      return result;
  }


// --------------- Printing to stream ---------------

/*
 * send each element of Vector to stream, seperated by a single space
 */
   template<typename Vec>
      requires is_vectorspace_base<Vec>::value
   std::ostream& operator<<( std::ostream& os, const Vec& v )
  {
      constexpr int N=Vec::N;

      for( int i=0; i<N-1; i++ )
     {
         os << v[i] << " ";
     }
      os << v[N-1];
      return os;
  }


