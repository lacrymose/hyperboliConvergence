
# ifndef ARRAY1D_H
# define ARRAY1D_H

# include <types.h>

# include <cassert>
# include <vector>
# include <algorithm>
# include <functional>

namespace Array
{

/*
 * Wrapper for std::vector providing element-wise arithmetic operations
 */
   template<typename T>
   struct Array1D
  {
      std::vector<T> data;

      Array1D( const int n=0 ){ data.resize(n); }

      void resize( const int n ){ data.resize(n); }
      size_t  size() const { return data.size(); }

      inline const T& operator[]( const int i ) const { return data[i]; }
      inline       T& operator[](       int i )       { return data[i]; }

      inline Array1D<T>& operator+=( const Array1D<T>& a );
      inline Array1D<T>& operator-=( const Array1D<T>& a );

      template<typename T2>
      inline Array1D<T>& operator+=( const Array1D<T2>& da );
      template<typename T2>
      inline Array1D<T>& operator-=( const Array1D<T2>& da );

      inline Array1D<T>& operator =(       Types::Real d );
      inline Array1D<T>& operator*=(       Types::Real d );
      inline Array1D<T>& operator/=(       Types::Real d );
  };

// Arithmetic operations

   template<typename T>
   inline Array1D<T> operator+( const Array1D<T>& a0, const Array1D<T>& a1 );
   template<typename T>
   inline Array1D<T> operator-( const Array1D<T>& a0, const Array1D<T>& a1 );

   template<typename T>
   inline Array1D<T> operator*( const Array1D<T>& a0,            Types::Real   d );
   template<typename T>
   inline Array1D<T> operator*(            Types::Real   d, const Array1D<T>& a1 );

   template<typename T>
   inline Array1D<T> operator/( const Array1D<T>& a0,            Types::Real   d );
}

# include <array1D/array1D.ipp>

# endif
