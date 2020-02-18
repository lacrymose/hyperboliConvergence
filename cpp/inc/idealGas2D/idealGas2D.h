# ifndef IDEAL_GAS_2D_H
# define IDEAL_GAS_2D_H

# include <ostream>

namespace IdealGas2D
{

/*
 * Vector of conserved variables: density, momentum (x,y) and total energy
 */
   struct ConservedVariables
  {
      float    var[4];

      inline ConservedVariables();

      inline float& operator[]( int i ){ return var[i]; }

      inline ConservedVariables  operator+(  ConservedVariables q );
      inline ConservedVariables  operator-(  ConservedVariables q );

      inline ConservedVariables& operator+=( ConservedVariables q );
      inline ConservedVariables& operator-=( ConservedVariables q );
      inline ConservedVariables& operator*=( float d );
      inline ConservedVariables& operator/=( float d );
      inline ConservedVariables& operator =( float d );

      inline friend std::ostream& operator<<( std::ostream& os, ConservedVariables q );
  };

/*
 * Vector of viscous variables: velocity(x,y), temperature and pressure
 */
   struct ViscousVariables
  {
      float    var[4];

      inline ViscousVariables();

      inline float& operator[]( int i ){ return var[i]; }

      inline ViscousVariables  operator+(  ViscousVariables q );
      inline ViscousVariables  operator-(  ViscousVariables q );

      inline ViscousVariables& operator+=( ViscousVariables q );
      inline ViscousVariables& operator-=( ViscousVariables q );
      inline ViscousVariables& operator*=( float d );
      inline ViscousVariables& operator/=( float d );
      inline ViscousVariables& operator =( float d );

      inline friend std::ostream& operator<<( std::ostream& os, ViscousVariables q );
  };

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

      inline State(){};
      inline State( Species& gas, ConservedVariables& qc );
      inline State( Species& gas, ViscousVariables&   qv );

      inline float& velocityX(){              return state[0]; }
      inline float& velocityY(){              return state[1]; }
      inline float& pressure(){               return state[2]; }
      inline float& density(){                return state[3]; }
      inline float& temperature(){            return state[4]; }
      inline float& specificTotalEnthalpy(){  return state[5]; }
      inline float& velocity2(){              return state[6]; }
      inline float& speedOfSound2(){          return state[7]; }
  };

/*
 * Nonlinear transformations for states
 */
   inline ConservedVariables conservedVariables( Species& gas, ConservedVariables& qc );
   inline ConservedVariables conservedVariables( Species& gas, ViscousVariables&   qv );
   inline ConservedVariables conservedVariables( Species& gas, State& s );

   inline ViscousVariables viscousVariables(     Species& gas, ConservedVariables& qc );
   inline ViscousVariables viscousVariables(     Species& gas, ViscousVariables&   qv );
   inline ViscousVariables viscousVariables(     Species& gas, State& s );

/*
 * Linear transformation of delta from conservative to viscous
 */
   inline ViscousVariables     dViscousVariables( Species& gas, State& state, ConservedVariables& dqc );

/*
 * Linear transformation of delta from viscous to conservative
 */
   inline ConservedVariables dConservedVariables( Species& gas, State& state,   ViscousVariables& dqv );

/*
 * Exact physical flux vector
 */
   inline void exactFlux(    Species& gas, float n[3], State& s,             ConservedVariables& f, float& lmax );

/*
 * Local Lax-Friedrichs flux
 */
   inline void laxFriedrichs( Species& gas, float n[3], State& sl, State& sr, ConservedVariables& f, float& lmax );

/*
 * AUSM+up flux for all speeds (M-S Liou 2006)
 */
   inline void ausm(         Species& gas, float n[3], State& sl, State& sr, ConservedVariables& f, float& lmax );
}

# include <idealGas2D/conservedVariables.ipp>
# include <idealGas2D/viscousVariables.ipp>
# include <idealGas2D/state.ipp>
# include <idealGas2D/transformations.ipp>
# include <idealGas2D/fluxes/ausm.ipp>
# include <idealGas2D/fluxes/exactFlux.ipp>
# include <idealGas2D/fluxes/laxFriedrichs.ipp>

# endif
