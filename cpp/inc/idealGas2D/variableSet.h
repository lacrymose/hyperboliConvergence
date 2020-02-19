# ifndef IDG2D_VARSET_H
# define IDG2D_VARSET_H

namespace IdealGas2D
{
/*
 * Vector of solution variables for an ideal gas, specified by template argument
 * Conservative <'c'>: ( density,  momentum,    total energy )
 * Viscous      <'v'>: ( velocity, temperature, pressure     )
 */
   template< char C >
   struct VariableSet
  {
      float var[4];

      inline VariableSet<C>();

      inline VariableSet<C>( const Species& gas, const VariableSet<C>& q0 );

      template< char D >
      inline VariableSet<C>( const Species& gas, const VariableSet<D>& q0 );

      inline VariableSet<C>( const Species& gas, const State& s0 );

      inline       float& operator[]( const int i )       { return var[i]; }
      inline const float& operator[]( const int i ) const { return var[i]; }

      inline VariableSet<C>& operator+=( const VariableSet<C>& q );
      inline VariableSet<C>& operator-=( const VariableSet<C>& q );

      inline VariableSet<C>& operator*=(       float d );
      inline VariableSet<C>& operator/=(       float d );
      inline VariableSet<C>& operator =(       float d );
  };

   template< char C >
   inline VariableSet<C> operator+( const VariableSet<C>& q0, const VariableSet<C>& q1 );
   template< char C >
   inline VariableSet<C> operator-( const VariableSet<C>& q0, const VariableSet<C>& q1 );

   template< char C >
   inline VariableSet<C> operator*( const VariableSet<C>& q0,                  float   d );
   template< char C >
   inline VariableSet<C> operator*(                float   d  , const VariableSet<C>& q1 );

   template< char C >
   inline VariableSet<C> operator/( const VariableSet<C>& q0,                  float   d );
}


# endif
