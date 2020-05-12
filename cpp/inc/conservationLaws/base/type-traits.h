
# pragma once

# include <conservationLaws/base/declarations.h>

# include <utils/type-traits.h>

# include <type_traits>

// ---------- dimension of types ----------

/*
 * returns the spatial dimension of a type, or -1 if it has no dimension parameter
 */
   template<typename T> struct dim_of;

// remove cvrefs
   template<typename T> struct dim_of<const T> : dim_of<T> {};
   template<typename T> struct dim_of<volatile T> : dim_of<T> {};
   template<typename T> struct dim_of<const volatile T> : dim_of<T> {};

   template<typename T> struct dim_of<T&> : dim_of<T> {};
   template<typename T> struct dim_of<T&&> : dim_of<T> {};


   // helpers
   template<typename T> constexpr int  dim_of_v = dim_of<T>::value;
   template<typename T> using dim_of_t = typename dim_of<T>::type;


/*
 * returns whether a type has a spatial dimension or not
 */
   template<typename T>
   struct has_dim : std::conditional_t< std::is_same_v< dim_of_t<T>,
                                                        std::integral_constant<int,-1> >,
                                        std::false_type,
                                        std::true_type> {};

   // helper
   template<typename T> constexpr bool has_dim_v = has_dim<T>::value;
   

/*
 * returns whether multiple types have the same spatial dimension
 */
   template<typename T0, typename... Ts>
   using has_same_dim = all< std::is_same< dim_of< remove_cvref_t<T0> >,
                                           dim_of< remove_cvref_t<Ts> > >... >;

   // helper
   template<typename T0, typename...Ts>
   constexpr bool has_same_dim_v = has_same_dim<T0,Ts...>::value;


/*
 * specialisations of dim_of to extract dimension from each type
 */

   // dim_of is -1 if template argument is not dimensioned
   template<typename T>
   struct dim_of : std::integral_constant<int,-1> {};

   template<LawType Law, int nDim>
   struct dim_of<State<Law,nDim>> : std::integral_constant<int,nDim> {};

   template<LawType Law, int nDim>
   struct dim_of<FluxResult<Law,nDim>> : std::integral_constant<int,nDim> {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct dim_of<VariableSet<Law,nDim,Basis>> : std::integral_constant<int,nDim> {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct dim_of<VariableDelta<Law,nDim,Basis>> : std::integral_constant<int,nDim> {};

/*
 * is dimension of type in {1,2,3}?
 */
   template<typename T>
      requires has_dim_v<T>
   struct has_valid_dim : std::bool_constant< ( (dim_of_v<T>) == 1 ) ||
                                              ( (dim_of_v<T>) == 2 ) ||
                                              ( (dim_of_v<T>) == 3 ) > {};

   // helper
   template<typename T>
   constexpr bool has_valid_dim_v = has_valid_dim<T>::value;


// ---------- law of types ----------

/*
 * returns the spatial law of a type, or -1 if it has no law parameter
 */
   template<typename T> struct law_of;

// remove cvrefs
   template<typename T> struct law_of<const T> : law_of<T> {};
   template<typename T> struct law_of<volatile T> : law_of<T> {};
   template<typename T> struct law_of<const volatile T> : law_of<T> {};

   template<typename T> struct law_of<T&> : law_of<T> {};
   template<typename T> struct law_of<T&&> : law_of<T> {};


   // helpers
   template<typename T> constexpr int  law_of_v = law_of<T>::value;
   template<typename T> using law_of_t = typename law_of<T>::type;


/*
 * returns whether a type has a spatial law or not
 */
   template<typename T>
   struct has_law : std::conditional_t< std::is_same_v< law_of_t<T>,
                                                        std::integral_constant<LawType,NoLaw> >,
                                        std::false_type,
                                        std::true_type> {};

   // helpers
   template<typename T> constexpr bool has_law_v = has_law<T>::value;
   template<typename T> using has_law_t = typename has_law<T>::type;
   

/*
 * returns whether multiple types have the same spatial law
 */
   template<typename T0, typename... Ts>
   using has_same_law = all< std::is_same< law_of< remove_cvref_t<T0> >,
                                           law_of< remove_cvref_t<Ts> > >... >;

   // helpers
   template<typename T0, typename...Ts>
   constexpr bool has_same_law_v = has_same_law<T0,Ts...>::value;


/*
 * specialisations of law_of to extract law from template types with typical parameter lists
 */

   // law_of is -1 if template argument is not lawed
   template<typename T>
   struct law_of : std::integral_constant<LawType,NoLaw> {};

   // covers all current use cases, and probably many future ones?
   template<template<LawType Law> T>
   struct law_of<T<Law>> : std::integral_constant<LawType,Law> {};

   template<template<LawType Law, int nDim> T>
   struct law_of<T<Law,nDim>> : std::integral_constant<LawType,Law> {};

   template<template<LawType Law, int nDim, BasisType<Law> Basis> T>
   struct law_of<T<Law,nDim,Basis>> : std::integral_constant<LawType,Law> {};


// ---------- return types of same law and spatial dimension ----------

   template<typename T> using    species_t =    Species< law_of_v<T> >;
   template<typename T> using      state_t =      State< law_of_v<T>, dim_of_v<T> >;
   template<typename T> using fluxresult_t = FluxResult< law_of_v<T>, dim_of_v<T> >;

