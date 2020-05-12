
# pragma once

// unique identifier for each conservation law
   enum struct LawType;


// number of scalar, vector or total variables required for each law
   template<LawType Law>
   constexpr int nScalarQuantities;

   template<LawType Law>
   constexpr int nVectorQuantities;

   template<LawType Law, int nDim>
   constexpr int nVar;


// typetraits to retrieve phase space basis enum for each law
   template<LawType Law> struct BasisTypeHelper;
   template<LawType Law> using BasisType = typename BasisTypeHelper<Law>::type;


// physical constants for each conservation law
   template<LawType Law > struct Species;

// unique point in phase space, agnostic to base
   template<LawType Law, int nDim> struct State;

// flux in conserved basis, and spectral radius
   template<LawType Law, int nDim> struct FluxResult;

// point in (affine) phase space with particular basis
   template<LawType Law, int nDim, BasisType<Law> Basis> struct VariableSet;

// increment in (affine) phase space with particular basis
   template<LawType Law, int nDim, BasisType<Law> Basis> struct VariableDelta;


