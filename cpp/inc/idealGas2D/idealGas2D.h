# ifndef IDEAL_GAS_2D_H
# define IDEAL_GAS_2D_H

# include <array1D/array1D.h>

# include <types.h>
# include <utils.h>

# include <array>
# include <iostream>
# include <iomanip>
# include <cassert>
# include <cmath>

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
// typedef VariableType<'p'>               Primitive;
// typedef VariableType<'g'>          Preconditioned;
// typedef VariableType<'w'>          Characteristic;

   typedef VariableSet<Conserved> ConservedVariables;
   typedef VariableSet<Viscous>     ViscousVariables;

   typedef VariableDelta<Conserved>   ConservedDelta;
   typedef VariableDelta<Viscous>       ViscousDelta;

/*
 * Species of ideal gas, including physical constants
 */
   struct Species
  {
      Types::Real gamma;
      Types::Real  Rgas;
      Types::Real    nu;
      Types::Real    pr;
      Types::Real  minf;

      void air();
  };

/*
 * Variety of data at a particular point in phase space, used for defining jacobians between different phase space bases or to hide particular choice of basis
 */
   struct State
  {
      std::array<Types::Real,8>  state;

      inline State();

   // create state using nonlinear transformation from VariableSet
      // delete prevents compiling for undefined VType
      template<typename VType>
      explicit inline State( const Species& gas, const VariableSet<VType>& q ) = delete;

      inline const Types::Real& velocityX()             const { return state[0]; }
      inline const Types::Real& velocityY()             const { return state[1]; }
      inline const Types::Real& pressure()              const { return state[2]; }
      inline const Types::Real& density()               const { return state[3]; }
      inline const Types::Real& temperature()           const { return state[4]; }
      inline const Types::Real& specificTotalEnthalpy() const { return state[5]; }
      inline const Types::Real& velocity2()             const { return state[6]; }
      inline const Types::Real& speedOfSound2()         const { return state[7]; }
  };

   inline bool operator==( const State& s0, const State& s1 );
   inline bool operator!=( const State& s0, const State& s1 );

/*
 * Vector of solution variables for an ideal gas, specified by template argument
 * Conservative <'c'>: ( density,  momentum,    total energy )
 * Viscous      <'v'>: ( velocity, temperature, pressure     )
 */
   template<typename VType>
   struct VariableSet
  {
   // variables
      std::array<Types::Real,4>  var;

   // default constructor
      inline VariableSet<VType>();

   // copy constructor
               inline VariableSet<VType>(                     const VariableSet<VType>&    q0 );

   // copy constructor to shadow nonlinear transformation from another set
      explicit inline VariableSet<VType>( const Species& gas, const VariableSet<VType>&    q0 );

   // convert dq -> q
      explicit inline VariableSet<VType>(                     const VariableDelta<VType>& dq0 );

   // nonlinear transformations from other variable sets
      // delete prevents compiling for undefined (VType,VType2) pairs
      template<typename VType2>
      explicit inline VariableSet<VType>( const Species& gas, const VariableSet<VType2>&   q0 ) = delete;
      explicit inline VariableSet<VType>( const Species& gas, const State&              state ) = delete;

   // accessors
      inline       Types::Real& operator[]( const int i )       { return var[i]; }
      inline const Types::Real& operator[]( const int i ) const { return var[i]; }

   // in-place arithmetic
      inline VariableSet<VType>& operator+=( const VariableDelta<VType>& dq0 );
      inline VariableSet<VType>& operator-=( const VariableDelta<VType>& dq0 );
      inline VariableSet<VType>& operator =(       Types::Real a );
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
      std::array<Types::Real,4>  var;

   // default constructor
      inline VariableDelta<VType>();

   // copy constructor
               inline VariableDelta<VType>(                     const VariableDelta<VType>& dq0 );

   // copy constructor to shadow linear transformation from another delta
      explicit inline VariableDelta<VType>( const Species& gas, const State& state, const VariableDelta<VType>& dq0 );

   // convert q -> dq
               inline VariableDelta<VType>(                     const VariableSet<VType>&    q0 );

   // linear transformations from other variable deltas
      template<typename VType2>
      explicit inline VariableDelta<VType>( const Species& gas, const State& state, const VariableDelta<VType2>& dq0 ) = delete;

   // accessors
      inline       Types::Real& operator[]( const int i )       { return var[i]; }
      inline const Types::Real& operator[]( const int i ) const { return var[i]; }

   // in-place arithmetic
      inline VariableDelta<VType>& operator+=( const VariableDelta<VType>& dq0 );
      inline VariableDelta<VType>& operator-=( const VariableDelta<VType>& dq0 );
      inline VariableDelta<VType>& operator*=(       Types::Real a );
      inline VariableDelta<VType>& operator/=(       Types::Real a );
      inline VariableDelta<VType>& operator =(       Types::Real a );
  };

   template<typename VType>
   inline VariableDelta<VType> abs( VariableDelta<VType> dq );

/*
 * Exact physical flux vector
 */
   inline void exactFlux(     const Species& gas, const Types::Real n[3], const State& state,               ConservedDelta& f, Types::Real& lmax );

/*
 * Local Rusanov flux
 */
   struct Rusanov
  {
      inline void operator()( const Species& gas, const Types::Real n[3], const State& sl, const State& sr, ConservedDelta& f, Types::Real& lmax ) const;
  };

/*
 * AUSM+up flux for all speeds (M-S Liou 2006)
 */
   struct Ausm
  {
      Types::Real alpha0;
      Types::Real  sigma;
      Types::Real   beta;
      Types::Real     Ku;
      Types::Real     Kp;

      Ausm();
      inline void operator()( const Species& gas, const Types::Real n[3], const State& sl, const State& sr, ConservedDelta& f, Types::Real& lmax ) const;
  };

/*
 * SLAU flux for all speeds (Shima & Kitamura 2011)
 */
   struct Slau
  {
      Types::Real alpha;

      Slau();
      inline void operator()( const Species& gas, const Types::Real n[3], const State& sl, const State& sr, ConservedDelta& f, Types::Real& lmax ) const;
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

# include <idealGas2D/fluxes/machSplittings.ipp>
# include <idealGas2D/fluxes/exactFlux.ipp>
# include <idealGas2D/fluxes/rusanov.ipp>
# include <idealGas2D/fluxes/ausm.ipp>
# include <idealGas2D/fluxes/slau.ipp>

# endif
