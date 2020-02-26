# ifndef UTILS_H
# define UTILS_H

# include <cmath>

namespace Utils
{
// return -1 if val<0, 0 if val==0, 1 if val>0
   template<typename T>
   inline int sign( const T val )
  {
      return ( T(0) < val ) - ( val < T(0) );
  }

// return 1 if val>0, else 0
   template<typename T>
   inline int greaterThan0( const T val )
  {
      return (sign(val)+1)/2;
  }

// return 1 if val<0, else 0
   template<typename T>
   inline int lessThan0( const T val )
  {
      return -(sign(val)-1)/2;
  }

// return 1 if val>off, else 0
   template<typename T>
   inline int greaterThan( const T val, const T offset )
  {
      return greaterThan0( val-offset );
  }

// return 1 if val<off, else 0
   template<typename T>
   inline int lessThan( const T val, const T offset )
  {
      return lessThan0( val-offset );
  }

// return val if val>0, else 0
   template<typename T>
   inline T positiveRamp( const T val )
  {
      return 0.5*( val + fabs( val ) );
  }

// return val if val<0, else 0
   template<typename T>
   inline T negativeRamp( const T val )
  {
      return 0.5*( val - fabs( val ) );
  }

// specializations for integral values
   template<>
   inline int positiveRamp( const int val )
  {
      return ( val + abs( val ) )/2;
  }

   template<>
   inline int negativeRamp( const int val )
  {
      return ( val - abs( val ) )/2;
  }
}

# endif
