# ifndef CONSERVATIONLAW_BASE_BASE_H
# define CONSERVATIONLAW_BASE_BASE_H

# include <types.h>

# include <array>

// forward declarations

   template<char C>
   struct EquationType;

   template<typename Law, char C>
   struct PhaseSpaceBasis;

   template<typename EqType, int NDIM>
   struct ConservationLaw;

   template<typename Law, typename Basis>
   struct VariableSet;

   template<typename Law, typename Basis>
   struct VariableDelta;

// typedefs

   typedef EquationType<'a'>        ScalarAdvection;
   typedef EquationType<'b'>        BurgersEquation;
   typedef EquationType<'e'>          EulerEquation;
   typedef EquationType<'s'>   ShallowWaterEquation;
   typedef EquationType<'i'>  IsothermalGasEquation;

/*
 * Blank type for unspecialises type aliases in class templates
 */
   struct NoType{};

/*
 * Equation Type describing how many scalar and vector quantities are required. Independent of number of physical dimensions
 *    e.g. the Euler Equations have one vector quantity (velocity/momentum etc) and two scalar quantities (thermodynamic states).
 *    nScalar Quantities must be specialised for each EquationType
 *    nVector Quantities must be specialised for each EquationType
 */
   template<char C>
   struct EquationType
  {
      static constexpr int nScalarQuantities;
      static constexpr int nVectorQuantities;
  };

/*
 * Label for a particular basis for the solution phase space, for example conserved or characteristic variables
 */
   template<typename Law, char C>
   PhaseSpaceBasis{};

/*
 * Static class holding values, types and functions for a generic hyperbolic conservation law
 *    Templated by an equation type, and the dimension of the physical space
 */
   template<typename EqType, int NDIM>
   struct ConservationLaw
  {
   // dimension of physical and phase spaces
      static constexpr int nDim = NDIM;
      static constexpr int nVar = EqType::nScalarQuantities + NDIM*EqType::nVectorQuantities;

   // these need to be specialised for each EqType
      using Species = NoType;
      using State   = NoType;

   // standard variable bases
      typedef PhaseSpaceBasis<ConservationLaw<EqType,nDim>,'c'>        ConservedBasis;
      typedef PhaseSpaceBasis<ConservationLaw<EqType,nDim>,'w'>   CharacteristicBasis;

   // conserved variable types
      typedef VariableSet<  ConservedBasis>            ConservedVariables;
      typedef VariableDelta<ConservedBasis>            ConservedDelta;

   // characteristic variable types
      typedef VariableSet<  CharacteristicBasis>  CharacteristicVariables;
      typedef VariableDelta<CharacteristicBasis>  CharacteristicDelta;

   // analytical expression for flux - needs to be specialised for each EqType
      static inline void exactFlux( const Species&                        species,
                                    const std::array<Types::Real,nDim+1>&       n,
                                    const State&                            state,
                                          ConservedDelta&                       f,
                                          Types::Real&                       lmax ) = delete;

      template<typename Basis>
      static inline void exactFlux( const Species&                                  species,
                                    const std::array<Types::Real,nDim+1>&                 n,
                                    const VariableSet<ConservationLaw<EqType,nDim>,Basis> q,
                                          ConservedDelta&                                 f,
                                          Types::Real&                                 lmax ) = delete;
  };

/*
 * Vector of solution variables for a hyperbolic conservation law (Law) in a particular variable set (Basis)
 *    Provides nonlinear transformations between VariableSets with the same Law but different Basis, and arithmetic operations with other VariableSets and VariableDeltas of the same Law and Basis
 */
   template<typename Law, typename Basis>
   struct VariableSet
  {
   // variables
      std::array<Types::Real,Law::nVar> var;

   // default constuctor
      inline VariableSet<Law,Basis>();

   // copy constructor
               inline VariableSet<Law,Basis>(                              const VariableSet<Law,Basis>&    q0 );

   // copy constructor to shadow nonlinear transformation from another set
      explicit inline VariableSet<Law,Basis>( const Law::Species& species, const VariableSet<Law,Basis>&    q0 );

   // convert dq -> q
      explicit inline VariableSet<Law,Basis>(                              const VariableDelta<Law,Basis>& dq0 );

   // nonlinear transformations from other variable sets
      // delete prevents compiling for unspecialised (Basis,Basis2) pairs
      template<typename Basis2>
      explicit inline VariableSet<Law,Basis>( const Law::Species& species, const VariableSet<Law,Basis2>&   q0 ) = delete;
      explicit inline VariableSet<Law,Basis>( const Law::Species& species, const Law::State&             state ) = delete;

   // accessors
      inline       Types::Real& operator[]( const int i )       { return var[i]; }
      inline const Types::Real& operator[]( const int i ) const { return var[i]; }

   // in-place arithmetic
      inline VariableSet<Law,Basis>& operator+=( const VariableDelta<Law,Basis>& dq0 );
      inline VariableSet<Law,Basis>& operator-=( const VariableDelta<Law,Basis>& dq0 );
      inline VariableSet<Law,Basis>& operator =( const Types::Real a );
  };

/*
 * Vector of solution variable deltas for a hyperbolic conservation law (Law) in a particular variable set (Basis)
 *    Provides linear transformations between VariableDeltas with the same Law but different Basis, and arithmetic operations with other VariableSets and VariableDeltas of the same Law and Basis
 */
   template<typename Law, typename Basis>
   struct VariableDelta
  {
   // variables
      std::array<Types::Real,Law::nVar> var;

   // default constructor
      inline VariableDelta<Law,Basis>();

   // copy constructor
               inline VariableDelta<Law,Basis>(                              const VariableDelta<Law,Basis>& dq0 );

   // copy constructor to shadow linear transformation from another delta
      explicit inline VariableDelta<Law,Basis>( const Law::Species& species, const Law::State& state, const VariableDelta<Law,Basis>& dq0 );

   // convert q -> dq
      explicit inline VariableDelta<Law,Basis>(                              const VariableSet<Law,Basis>& dq0 );

   // linear transformations from other variable deltas
      // delete prevents compiling for unspecialised (Basis,Basis2) pairs
      template<typename Basis2>
      explicit inline VariableDelta<Law,Basis>( const Law::Species& species, const Law::State& state, const VariableDelta<Law,Basis2>&   q0 ) = delete;

   // accessors
      inline       Types::Real& operator[]( const int i )       { return var[i]; }
      inline const Types::Real& operator[]( const int i ) const { return var[i]; }

   // in-place arithmetic
      inline VariableDelta<Law,Basis>& operator+=( const VariableDelta<Law,Basis>& dq0 );
      inline VariableDelta<Law,Basis>& operator-=( const VariableDelta<Law,Basis>& dq0 );
      inline VariableDelta<Law,Basis>& operator*=( const Types::Real a );
      inline VariableDelta<Law,Basis>& operator/=( const Types::Real a );
      inline VariableDelta<Law,Basis>& operator =( const Types::Real a );
  };

/*
 * Interface flux using a central average
 */
   template<typename Law>
   struct CentralFlux
  {
      inline void operator()( const Law::Species&                    species,
                              const std::array<Types::Real,Law::nDim+1>&   n,
                              const Law::State&                           sl,
                              const Law::State&                           sr,
                                    Law::ConservedDelta&                   f,
                                    Types::Real&                        lmax );

      template<typename Basis>
      inline void operator()( const Law::Species&                    species,
                              const std::array<Types::Real,Law::nDim+1>&   n,
                              const VariableSet<Law,Basis>&               ql,
                              const VariableSet<Law,Basis>&               qr,
                                    Law::ConservedDelta&                   f,
                                    Types::Real&                        lmax );
  };

/*
 * Interface using a Rusanov flux with scalar dissipation
 */
   template<typename Law>
   struct RusanovFlux
  {
      inline void operator()( const Law::Species&                    species,
                              const std::array<Types::Real,Law::nDim+1>&   n,
                              const Law::State&                           sl,
                              const Law::State&                           sr,
                                    Law::ConservedDelta&                   f,
                                    Types::Real&                        lmax );

      template<typename Basis>
      inline void operator()( const Law::Species&                    species,
                              const std::array<Types::Real,Law::nDim+1>&   n,
                              const VariableSet<Law,Basis>&               ql,
                              const VariableSet<Law,Basis>&               qr,
                                    Law::ConservedDelta&                   f,
                                    Types::Real&                        lmax );
  };

# include <conservationLaws/base/fluxes/centralFlux.ipp>
# include <conservationLaws/base/fluxes/RusanovFlux.ipp>

# endif
