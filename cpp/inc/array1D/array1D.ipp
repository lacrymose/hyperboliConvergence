
   template<typename T>
   inline Array1D<T>& Array1D<T>::operator =( float d)
  {
      for( T& t : data ){ t=d; }
      return *this;
  }

   template<typename T>
   inline Array1D<T>& Array1D<T>::operator*=( float d)
  {
      for( T& t : data ){ t*=d; }
      return *this;
  }

   template<typename T>
   inline Array1D<T>& Array1D<T>::operator/=( float d)
  {
      for( T& t : data ){ t/=d; }
      return *this;
  }

   template<typename T>
   inline Array1D<T>& Array1D<T>::operator+=( const Array1D<T>& a0 )
  {
      assert( size() == a0.size() );

      std::transform(    data.begin(), data.end(),
                      a0.data.begin(),
                         data.begin(),
                         std::plus<T>() );

      return *this;
  }

   template<typename T>
   inline Array1D<T>& Array1D<T>::operator-=( const Array1D<T>& a0 )
  {
      assert( size() == a0.size() );

      std::transform(    data.begin(), data.end(),
                      a0.data.begin(),
                         data.begin(),
                         std::minus<T>() );

      return *this;
  }

   template<typename T>
   inline Array1D<T> operator+( const Array1D<T>& a0, const Array1D<T>& a1 )
  {
      int n0=a0.size();
      int n1=a1.size();
      assert( n0 == n1 );

      Array1D<T>  a2( n0 );

      std::transform( a0.data.begin(), a0.data.end(),
                      a1.data.begin(),
                      a2.data.begin(),
                         std::plus<T>() );
      return a2;
  }

   template<typename T>
   inline Array1D<T> operator-( const Array1D<T>& a0, const Array1D<T>& a1 )
  {
      int n0=a0.size();
      int n1=a1.size();
      assert( n0 == n1 );

      Array1D<T>  a2( n0 );

      std::transform( a0.data.begin(), a0.data.end(),
                      a1.data.begin(),
                      a2.data.begin(),
                         std::minus<T>() );
      return a2;
  }

   template<typename T>
   inline Array1D<T> operator*( const Array1D<T>& a0, float d )
  {
      int n=a0.size();
      Array1D<T>  a1( n );
      for( int i=0; i<n; i++ )
     {
         a1[i] = a0[i]*d;
     }
      return a1;
  }

   template<typename T>
   inline Array1D<T> operator*( float d, const Array1D<T>& a0 )
  {
      int n=a0.size();
      Array1D<T>  a1( n );
      for( int i=0; i<n; i++ )
     {
         a1[i] = d*a0[i];
     }
      return a1;
  }

   template<typename T>
   inline Array1D<T> operator/( const Array1D<T>& a0, float d )
  {
      int   n=a0.size();
      float     d1=1./d;

      Array1D<T>  a1( n );
      for( int i=0; i<n; i++ )
     {
         a1[i] = a0[i]*d1;
     }
      return a1;
  }

