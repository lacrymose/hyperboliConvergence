
# pragma once

# include <conservationLaws/base/type-traits/dim.h>
# include <conservationLaws/base/type-traits/law.h>
# include <conservationLaws/base/type-traits/base.h>

# include <conservationLaws/base/declarations.h>

# include <utils/type-traits.h>

# include <type_traits>

// ---------- return types of same law and spatial dimension ----------

   template<typename T> using    species_t =    Species< law_of_v<T> >;
   template<typename T> using      state_t =      State< law_of_v<T>, dim_of_v<T> >;
   template<typename T> using fluxresult_t = FluxResult< law_of_v<T>, dim_of_v<T> >;


// ---------- is type an instantiation of a particular type template? ----------

// Species ----------

   template<typename T>
   struct is_Species : std::false_type {};

   template<LawType Law>
   struct is_Species<Species<Law>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_Species_v = is_Species<T>::value;


// State ----------

   template<typename T>
   struct is_State : std::false_type {};

   template<LawType Law, int nDim>
   struct is_State<State<Law,nDim>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_State_v = is_State<T>::value;


// FluxResult ----------

   template<typename T>
   struct is_FluxResult : std::false_type {};

   template<LawType Law, int nDim>
   struct is_FluxResult<FluxResult<Law,nDim>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_FluxResult_v = is_FluxResult<T>::value;


// VariableSet ----------

   template<typename T>
   struct is_VariableSet : std::false_type {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct is_VariableSet<VariableSet<Law,nDim,Basis>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_VariableSet_v = is_VariableSet<T>::value;


// VariableDelta ----------

   template<typename T>
   struct is_VariableDelta : std::false_type {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct is_VariableDelta<VariableDelta<Law,nDim,Basis>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_VariableDelta_v = is_VariableDelta<T>::value;


// VariableSet with specific Law and Basis ----------

   template<typename T, LawType Law, BasisType<Law> Basis>
   struct is_specialised_VarSet : std::false_type {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct is_specialised_VarSet<VariableSet<Law,nDim,Basis>,Law,Basis> : std::true_type{};

   // helper
   template<typename T, LawType Law, BasisType<Law> Basis>
   constexpr bool is_specialised_VarSet_v = is_specialised_VarSet<T,Law,Basis>::value;


// VariableDelta with specific Law and Basis ----------

   template<typename T, LawType Law, BasisType<Law> Basis>
   struct is_specialised_VarDelta : std::false_type {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct is_specialised_VarDelta<VariableDelta<Law,nDim,Basis>,Law,Basis> : std::true_type{};

   // helper
   template<typename T, LawType Law, BasisType<Law> Basis>
   constexpr bool is_specialised_VarDelta_v = is_specialised_VarDelta<T,Law,Basis>::value;

