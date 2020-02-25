# ifndef IDEAL_GAS_2D_H
# define IDEAL_GAS_2D_H

# include <array1D/array1D.h>

# include <iostream>
# include <iomanip>
# include <assert.h>
# include <math.h>

namespace IdealGas2D
{
   template<char C>
   struct VariableType{};

// forward declarations
   template<typename VType>
   struct VariableSet;

   template<typename VType>
   struct VariableDelta;

   struct Species;
   struct State;

// typedefs
   typedef VariableType<'c'>               Conserved;
   typedef VariableType<'v'>                 Viscous;

   typedef VariableSet<Conserved> ConservedVariables;
   typedef VariableSet<Viscous>     ViscousVariables;

   typedef VariableDelta<Conserved>   ConservedDelta;
   typedef VariableDelta<Viscous>       ViscousDelta;

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

   // create state using nonlinear transformation from VariableSet
      // delete prevents compiling for undefined VType
      template<typename VType>
      explicit inline State( const Species& gas, const VariableSet<VType>& q ) = delete;

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
   template<typename VType>
   struct VariableSet
  {
   // variables
      float var[4];

   // default constructor
      inline VariableSet<VType>();

   // copy constructors
               inline VariableSet<VType>(                     const VariableSet<VType>&    q0 );
      explicit inline VariableSet<VType>( const Species& gas, const VariableSet<VType>&    q0 );

   // convert dq -> q
      explicit inline VariableSet<VType>(                     const VariableDelta<VType>& dq0 );

   // nonlinear transformations from other variable sets
      // delete prevents compiling for undefined (VType,VType2) pairs
      template<typename VType2>
      explicit inline VariableSet<VType>( const Species& gas, const VariableSet<VType2>&   q0 ) = delete;
      explicit inline VariableSet<VType>( const Species& gas, const State&              state ) = delete;

   // accessors
      inline       float& operator[]( const int i )       { return var[i]; }
      inline const float& operator[]( const int i ) const { return var[i]; }

   // in-place arithmetic
      inline VariableSet<VType>& operator+=( const VariableDelta<VType>& dq0 );
      inline VariableSet<VType>& operator-=( const VariableDelta<VType>& dq0 );
      inline VariableSet<VType>& operator =(       float a );
  };

/*
 * Vector of solution variable deltas for an ideal gas, specified by template argument
 * Conservative <'c'>: ( density,  momentum,    total energy )
 * Viscous      <'v'>: ( velocity, temperature, pressure     )
 */
   template<typename VType>
   struct VariableDelta
  {
   // variables
      float var[4];

   // default constructor
      inline VariableDelta<VType>();

   // copy constructors
               inline VariableDelta<VType>(                     const VariableDelta<VType>& dq0 );
      explicit inline VariableDelta<VType>( const Species& gas, const VariableDelta<VType>& dq0 );

   // convert q -> dq
               inline VariableDelta<VType>(                     const VariableSet<VType>&    q0 );

   // linear transformations from other variable deltas
      template<typename VType2>
      explicit inline VariableDelta<VType>( const Species& gas, const State& state, const VariableDelta<VType2>& dq0 ) = delete;

   // accessors
      inline       float& operator[]( const int i )       { return var[i]; }
      inline const float& operator[]( const int i ) const { return var[i]; }

   // in-place arithmetic
      inline VariableDelta<VType>& operator+=( const VariableDelta<VType>& dq0 );
      inline VariableDelta<VType>& operator-=( const VariableDelta<VType>& dq0 );
      inline VariableDelta<VType>& operator*=(       float a );
      inline VariableDelta<VType>& operator/=(       float a );
      inline VariableDelta<VType>& operator =(       float a );
  };

/*
 * Exact physical flux vector
 */
   inline void exactFlux(     const Species& gas, const float n[3], const State& state,               ConservedDelta& f, float& lmax );

/*
 * Local Lax-Friedrichs flux
 */
   struct LaxFriedrichs
  {
      inline void operator()( const Species& gas, const float n[3], const State& sl, const State& sr, ConservedDelta& f, float& lmax ) const;
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
      inline void operator()( const Species& gas, const float n[3], const State& sl, const State& sr, ConservedDelta& f, float& lmax ) const;
  };
}

# include <idealGas2D/state.ipp>

# include <idealGas2D/variables/variableSet.ipp>
# include <idealGas2D/variables/variableDelta.ipp>
# include <idealGas2D/variables/variableArithmetic.ipp>

# include <idealGas2D/variables/viscousSet.ipp>
# include <idealGas2D/variables/viscousDelta.ipp>

# include <idealGas2D/variables/conservedSet.ipp>
# include <idealGas2D/variables/conservedDelta.ipp>

# include <idealGas2D/fluxes/exactFlux.ipp>
# include <idealGas2D/fluxes/ausm.ipp>
# include <idealGas2D/fluxes/laxFriedrichs.ipp>

# endif
