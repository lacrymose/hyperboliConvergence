# ifndef IDEAL_GAS_2D_H
# define IDEAL_GAS_2D_H

# include <iostream>
# include <assert.h>
# include <math.h>


namespace IdealGas2D
{
// forward declarations
   struct Species;
   struct State;

   template<char C>
   struct VariableSet;

// typedefs
   typedef VariableSet<'c'> ConservedVariables;
   typedef VariableSet<'v'> ViscousVariables;

/*
 * Species of ideal gas, including physical constants
 */
   struct Species
  {
      float gamma;
      float  Rgas;
      float    nu;
      float    pr;
      float  minf;

      void air();
  };

/*
 * Particular point in phase space, used for defining jacobians between different phase space bases
 */
   struct State
  {
      float    state[8];

      inline State();

      template<char C>
      explicit inline State( const Species& gas, const VariableSet<C>& q );

      inline const float& velocityX()             const { return state[0]; }
      inline const float& velocityY()             const { return state[1]; }
      inline const float& pressure()              const { return state[2]; }
      inline const float& density()               const { return state[3]; }
      inline const float& temperature()           const { return state[4]; }
      inline const float& specificTotalEnthalpy() const { return state[5]; }
      inline const float& velocity2()             const { return state[6]; }
      inline const float& speedOfSound2()         const { return state[7]; }
  };

/*
 * Vector of solution variables for an ideal gas, specified by template argument
 * Conservative <'c'>: ( density,  momentum,    total energy )
 * Viscous      <'v'>: ( velocity, temperature, pressure     )
 */
   template<char C>
   struct VariableSet
  {
      float var[4];

      inline VariableSet<C>();

               inline VariableSet<C>(                     const VariableSet<C>& q0 );

      explicit inline VariableSet<C>( const Species& gas, const VariableSet<C>& q0 );

      template<char D>
      explicit inline VariableSet<C>( const Species& gas, const VariableSet<D>& q0 );

      explicit inline VariableSet<C>( const Species& gas, const State& s0 );

      inline       float& operator[]( const int i )       { return var[i]; }
      inline const float& operator[]( const int i ) const { return var[i]; }

      inline VariableSet<C>& operator+=( const VariableSet<C>& q );
      inline VariableSet<C>& operator-=( const VariableSet<C>& q );

      inline VariableSet<C>& operator*=(       float d );
      inline VariableSet<C>& operator/=(       float d );
      inline VariableSet<C>& operator =(       float d );
  };

// friend functions for VariableSet
   template<char C>
   inline VariableSet<C> operator+( const VariableSet<C>& q0, const VariableSet<C>& q1 );
   template<char C>
   inline VariableSet<C> operator-( const VariableSet<C>& q0, const VariableSet<C>& q1 );

   template<char C>
   inline VariableSet<C> operator*( const VariableSet<C>& q0,                  float   d );
   template<char C>
   inline VariableSet<C> operator*(                float   d  , const VariableSet<C>& q1 );

   template<char C>
   inline VariableSet<C> operator/( const VariableSet<C>& q0,                  float   d );

/*
 * Exact physical flux vector
 */
   inline void exactFlux(     const Species& gas, const float n[3], const State& s,                   VariableSet<'c'>& f, float& lmax );

/*
 * Local Lax-Friedrichs flux
 */
   struct LaxFriedrichs
  {
    inline void operator()( const Species& gas, const float n[3], const State& sl, const State& sr, VariableSet<'c'>& f, float& lmax ) const;
  };

/*
 * AUSM+up flux for all speeds (M-S Liou 2006)
 */
   struct Ausm
  {
      float alpha0;
      float  sigma;
      float   beta;
      float     Ku;
      float     Kp;

      Ausm();
      inline void operator()( const Species& gas, const float n[3], const State& sl, const State& sr, VariableSet<'c'>& f, float& lmax ) const;
  };

}

# include <idealGas2D/state.ipp>
# include <idealGas2D/variables/variableSet.ipp>
# include <idealGas2D/variables/conservedVariables.ipp>
# include <idealGas2D/variables/viscousVariables.ipp>
# include <idealGas2D/fluxes/exactFlux.ipp>
# include <idealGas2D/fluxes/ausm.ipp>
# include <idealGas2D/fluxes/laxFriedrichs.ipp>

# endif
