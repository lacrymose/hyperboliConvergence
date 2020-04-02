# ifndef CONSERVATIONLAW_BASE_BASE_H
# define CONSERVATIONLAW_BASE_BASE_H

# include <types.h>

# include <geometry.h>

# include <array>

// forward declarations

   template<char LawLabel>
   struct ConservationLaw;

   template<typename Law>
   struct Species;

   template<typename Law, int nDim>
   struct State;

   template<typename Law, char BasisLabel, int nDim>
   struct PhaseSpaceBasis;

   template<typename Law, typename Basis>
   struct VariableSet;

   template<typename Law, typename Basis>
   struct VariableDelta;

   template<typename Law, int nDim>
   struct FluxResult;

// typedefs

   static constexpr char            ScalarAdvectionLabel = 'a';
   static constexpr char                    BurgersLabel = 'b';
   static constexpr char  ArtificialCompressibilityLabel = 'c';
   static constexpr char                      EulerLabel = 'e';
   static constexpr char              IsothermalGasLabel = 'i';
   static constexpr char               ShallowWaterLabel = 's';
   static constexpr char                TrafficFlowLabel = 't';

   static constexpr char             ConservedBasisLabel = 'c';
   static constexpr char        CharacteristicBasisLabel = 'w';

/*
 * Static class holding values, types and functions for a generic hyperbolic conservation law
 *       LawLabel must be the label associated with an implemented type of Law
 */
   template<char LawLabel>
   struct ConservationLaw
  {
   // number and type of conserved quantities
   //    must be specialised for each LawLabel
      static constexpr int  nScalarQuantities;
      static constexpr int  nVectorQuantities;

   // total number of variables (dimension of phase space)
      template<int nDim>
      static constexpr int  nVar = nScalarQuantities + nDim*nVectorQuantities;

   // species of conserved quantities - contains physical properties
      using ::Species<ConservationLaw>  = Species;

   // Basis-independent state
      template<nDim>
      using ::State<ConservationLaw,nDim> = State;

   // specialised phase space bases
      template<char BasisLabel, int nDim>
      using ::PhaseSpaceBasis<ConservationLaw,BasisLabel,nDim> = PhaseSpaceBasis;

   // specialised variable vectors
      template<typename Basis> using ::VariableSet<  ConservationLaw,Basis> = VariableSet;
      template<typename Basis> using ::VariableDelta<ConservationLaw,Basis> = VariableDelta;

   // standard variable bases
      // conserved variables
      template<int nDim> using PhaseSpaceBasis<'c',nDim>           = ConservedBasis;
      template<int nDim> using VariableSet<  ConservedBasis<nDim>> = ConservedVariables;
      template<int nDim> using VariableDelta<ConservedBasis<nDim>> = ConservedDelta;

      // characteristic variables
      template<int nDim> using PhaseSpaceBasis<'w',nDim>                = CharacteristicBasis;
      template<int nDim> using VariableDelta<CharacteristicBasis<nDim>> = CharacteristicDelta;

   // specialised return values of flux functions
      template<int nDim> using ::FluxResult<ConservationLaw,nDim> = FluxResult;

   // specialised standard flux functors
      using ::CentralFLux<ConservationLaw> = CentralFlux;
      using ::RusanovFLux<ConservationLaw> = RusanovFlux;

   // analytical expression for flux
   //    must be specialised for each LawLabel
      template<int nDim>
      static FluxResult<nDim> exactFlux( const Species              species,
                                         const Geometry::Quad<nDim>    face,
                                         const State<nDim>            state ) = delete;

   // analytical expression for flux
   //    this defaults to a wrapper calling exactFlux with a State argument but can be specialised for faster performance
      template<typename Basis>
      static FluxResult<Basis::nDim> exactFlux( const Species                     species,
                                                const Geometry::Quad<Basis::nDim>    face,
                                                const VariableSet<Basis>                q );
  };

/*
 * Physical constants
 *    Includes properties such as ratio of specific heats for ideal gas, gravity for shallow water equations etc
 *       Law must be an instantiated ConservationLaw template class
 */
   template<typename Law>
   struct Species{};

/*
 * Point in phase space independent of any phase space basis.
 *    Used to hide particular phase space basis, and save common computations
 *       Law must be an instantiated ConservationLaw template class
 */
   template<typename Law, int nDim>
   struct State
  {
   // state values
   //    number of elements must be specialised for each Law
      std::array<Types::Real,Law::nVar<nDim>> state;

   // default constructor
      State();

   // copy constructor
      State( State s0 );

   // nonlinear transformation from a variable set
   //    must be specialised for each Basis
      template<typename Basis>
      explicit State( const Species<Law>           species,
                      const VariableSet<Law,Basis>       q ) = delete;

   // accessors with verbose names for retrieving elements of state
   //    must be specialised for each Law

   // Types::Real some_state_element() const { return state[some_index]; }
  };

/*
 * Basis for the solution phase space, for example conserved or characteristic variables
 *       Law must be an instantiated ConservationLaw template class
 *       BasisLabel determines the particular basis
 *       NDIM is the dimension of the physical space: {1,2,3}
 */
   template<typename Law, char BasisLabel, int NDIM>
   struct PhaseSpaceBasis
  {
      static constexpr int nDim=NDIM;
      static constexpr int nVar=Law::nVar<nDim>;
  };

/*
 * Vector of solution variables for a hyperbolic conservation law (Law) in a particular variable set (Basis)
 *    Provides nonlinear transformations between VariableSets with the same Law but different Basis, and arithmetic operations with other VariableSets and VariableDeltas of the same Law and Basis
 *       Law   must be an instantiated ConservationLaw template class
 *       Basis must be an instantiated PhaseSpaceBasis template class
 */
   template<typename Law, typename Basis>
   struct VariableSet
  {
   // variables
      std::array<Types::Real,Basis::nVar> var;

   // default constructor
      VariableSet();

   // copy constructor
      VariableSet( VariableSet q0 );

   // convert dq -> q
      explicit VariableSet( const VariableDelta<Law,Basis> dq );

   // nonlinear transformation from another basis
   //    defaults to converting Basis2 -> State -> Basis but can be specialized for faster performance
      template<typename Basis2>
      explicit VariableSet( const Species<Law>            species,
                            const VariableSet<Law,Basis2>      q0 );

   // nonlinear transformation from a state
   //    must be specialised for each {Law,Basis} pair
      explicit VariableSet( const Species<Law>           species,
                            const State<Law,Basis::nDim>  state ) = delete;

   // accessors
            Types::Real& operator[]( const int i );
      const Types::Real& operator[]( const int i ) const;

   // in-place arithmetic
      VariableSet& operator+=( const VariableDelta<Law,Basis> dq0 );
      VariableSet& operator-=( const VariableDelta<Law,Basis> dq0 );
      VariableSet& operator =( const Types::Real a );
  };

/*
 * Vector of solution variable deltas for a hyperbolic conservation law (Law) in a particular variable set (Basis)
 *    Provides linear transformations between VariableDeltas with the same Law but different Basis, and arithmetic operations with other VariableSets and VariableDeltas of the same Law and Basis
 *       Law   must be an instantiated ConservationLaw template class
 *       Basis must be an instantiated PhaseSpaceBasis template class
 */
   template<typename Law, typename Basis>
   struct VariableDelta
  {
   // variables
      std::array<Types::Real,Basis::nVar> var;

   // default constructor
      VariableDelta();

   // copy constructor
      VariableDelta( VariableDelta q0 );

   // convert q -> dq
      explicit VariableDelta( const VariableSet<Law,Basis> q );

   // linear transformation from another basis
   //    this needs to be specialised for each {Law,Basis,Basis2} set
      template<typename Basis2>
      explicit VariableDelta( const Species<Law>              species,
                              const State<Law,Basis::nDim>      state,
                              const VariableDelta<Law,Basis2>     dq0 ) = delete;

   // accessors
            Types::Real& operator[]( const int i );
      const Types::Real& operator[]( const int i ) const;

   // in-place arithmetic
      VariableDelta& operator+=( const VariableDelta dq0 );
      VariableDelta& operator-=( const VariableDelta dq0 );
      VariableDelta& operator*=( const Types::Real a );
      VariableDelta& operator/=( const Types::Real a );
      VariableDelta& operator =( const Types::Real a );
  };

/*
 * Flux returns a delta of conserved values, and a spectral radius
 *       Law  must be an instantiated ConservationLaw template class
 *       NDIM is the dimension of the physical space: {1,2,3}
 */
   template<typename Law, int nDim>
   struct FluxResult
  {
      typedef typename Law::ConservedDelta<nDim> FluxType;

      FluxType      flux;
      Types::Real lambda;
  };

/*
 * face flux interface
 *       Law  must be an instantiated ConservationLaw template class
 *       Flux must have a method named flux with the call signature Species, Quad, State, State
 */
   template<typename Law>
   struct FluxInterface
  {
      template<typename Basis>
      FluxResult<Law,Basis::nDim> operator()( const Species<Law>                species,
                                              const Geometry::Quad<Basis::nDim>    face,
                                              const VariableSet<Basis>               ql,
                                              const VariableSet<Basis>               qr );

      template<int nDim>
      FluxResult<Law,nDim> operator()( const Species<Law>         species,
                                       const Geometry::Quad<nDim>    face,
                                       const State<Law,nDim>           sl,
                                       const State<Law,nDim>           sr );
  };

/*
 * Interface flux using a central average
 *       Law  must be an instantiated ConservationLaw template class
 */
   template<typename Law>
   struct CentralFlux : FluxInterface<Law>
  {
      template<int nDim>
      static FluxResult<Law,nDim> flux( const Species<Law>         species,
                                        const Geometry::Quad<nDim>    face,
                                        const State<Law,nDim>           sl,
                                        const State<Law,nDim>           sr );
  };

/*
 * Interface using a Rusanov flux with scalar dissipation
 *       Law  must be an instantiated ConservationLaw template class
 */
   template<typename Law>
   struct RusanovFlux : FluxInterface<Law>
  {
      template<int nDim>
      static FluxResult<Law,nDim> flux( const Species<Law>         species,
                                        const Geometry::Quad<nDim>    face,
                                        const State<Law,nDim>           sl,
                                        const State<Law,nDim>           sr );
  };

# include <conservationLaws/base/fluxes/exactFlux.ipp>
# include <conservationLaws/base/fluxes/centralFlux.ipp>
# include <conservationLaws/base/fluxes/RusanovFlux.ipp>

# endif
