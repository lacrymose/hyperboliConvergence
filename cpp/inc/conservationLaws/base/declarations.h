
# pragma once

# include <utils/concepts.h>

// unique identifier for each conservation law
   enum struct LawType { NoLaw,
                         ArtificialCompressibility,
                         Burgers,
                         Euler,
                         IsothermalGas,
                         Maxwell,
                         ScalarAdvection,
                         ShallowWater,
                         TrafficFlow };


// number of scalar, vector or total variables required for each law
   template<LawType Law>
   constexpr int nScalarQuantities{0};

   template<LawType Law>
   constexpr int nVectorQuantities{0};

   template<LawType Law, int nDim>
   constexpr int nVar = nScalarQuantities<Law> + nDim*nVectorQuantities<Law>;


// typetraits to retrieve phase space basis enum for each law
   template<LawType Law> struct BasisTypeHelper;
   template<LawType Law> using  BasisType = typename BasisTypeHelper<Law>::type;


// physical constants for each conservation law
   template<LawType Law, floating_point Real> struct Species;

// unique point in phase space, agnostic to base
   template<LawType Law, int nDim, floating_point Real> struct State;

// flux in conserved basis, and spectral radius
   template<LawType Law, int nDim, floating_point Real> struct FluxResult;

// point in (affine) phase space with particular basis
   template<LawType Law, int nDim, BasisType<Law> Basis, floating_point Real> struct VariableSet;

// increment in (affine) phase space with particular basis
   template<LawType Law, int nDim, BasisType<Law> Basis, floating_point Real> struct VariableDelta;


// transformations between variable sets and states
// template<typename DstT, typename SrcT, typename SpeciesT>
// DstT set2Set( const SpeciesT&, const SrcT& )=delete;

   template<typename DstT, typename StateT, typename SrcT, typename SpeciesT>
   DstT delta2Delta( const SpeciesT&, const StateT&, const SrcT& )=delete;

   template<typename DstT, typename SrcT, typename SpeciesT>
   auto set2State( const SpeciesT&, const SrcT& )=delete;

   template<typename DstT, typename SrcT, typename SpeciesT>
   DstT state2Set( const SpeciesT&, const SrcT& )=delete;

