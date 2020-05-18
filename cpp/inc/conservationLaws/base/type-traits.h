
# pragma once

# include <conservationLaws/base/declarations.h>

# include <utils/type-traits.h>

# include <type_traits>

// ---------- dimension of types ----------

// turns a dimension value into a type so it can be used in has_dim etc
   template<int nDim>
   using dim_constant = std::integral_constant<int,nDim>;


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
   struct has_dim : std::conditional_t< dim_of_v<T> == -1,
                                        std::false_type,
                                        std::true_type> {};

   // helper
   template<typename T> constexpr bool has_dim_v = has_dim<T>::value;
   

/*
 * returns whether multiple types have the same spatial dimension
 */
   template<typename T0, typename... Ts>
   using has_same_dim = all<std::is_same<dim_of_t<T0>,
                                         dim_of_t<Ts>>...>;

   // helper
   template<typename T0, typename...Ts>
   constexpr bool has_same_dim_v = has_same_dim<T0,Ts...>::value;


/*
 * specialisations of law_of to extract law from template types with typical parameter lists
 */

   // dim_of is -1 if template argument is not dimensioned
   template<typename T>
   struct dim_of : dim_constant<-1> {};

   template<int nDim>
   struct dim_of<dim_constant<nDim>> : dim_constant<nDim> {};

   // covers all current use cases, and probably many future ones?
   template<LawType Law, int nDim, template<LawType, int> typename T>
   struct dim_of<T<Law,nDim>> : dim_constant<nDim> {};

   template<LawType Law, int nDim, BasisType<Law> Basis, template<LawType, int, BasisType<Law>> typename T>
   struct dim_of<T<Law,nDim,Basis>> : dim_constant<nDim> {};

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

// turns a LawType value into a type so it can be used in has_law etc
   template<LawType Law>
   using law_constant = std::integral_constant<LawType,Law>;

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
   template<typename T> constexpr LawType law_of_v = law_of<T>::value;
   template<typename T> using law_of_t = typename law_of<T>::type;


/*
 * returns whether a type has a spatial law or not
 */
   template<typename T>
   struct has_law : std::conditional_t< std::is_same_v< law_of_t<T>,
                                                        law_constant<LawType::NoLaw> >,
                                        std::false_type,
                                        std::true_type> {};

// template<typename T>
// struct has_law : std::conditional_t< (law_of_v<T>==LawType::NoLaw),
//                                      std::false_type,
//                                      std::true_type> {};

   // helpers
   template<typename T> constexpr bool has_law_v = has_law<T>::value;
// template<typename T> using has_law_t = typename has_law<T>::type;
   

/*
 * returns whether multiple types have the same spatial law
 */
   template<typename T0, typename... Ts>
   using has_same_law = all< std::is_same<law_of_t<T0>,
                                          law_of_t<Ts>>... >;

   // helpers
   template<typename T0, typename...Ts>
   constexpr bool has_same_law_v = has_same_law<T0,Ts...>::value;


/*
 * specialisations of law_of to extract law from template types with typical parameter lists
 */

   // law_of is -1 if template argument is not lawed
   template<typename T>
   struct law_of : law_constant<LawType::NoLaw> {};

   template<LawType Law>
   struct law_of<law_constant<Law>> : law_constant<Law> {};

   // covers all current use cases, and probably many future ones?
   template<LawType Law, template<LawType> typename T>
   struct law_of<T<Law>> : law_constant<Law> {};

   template<LawType Law, int nDim, template<LawType, int> typename T>
   struct law_of<T<Law,nDim>> : law_constant<Law> {};

   template<LawType Law, int nDim, BasisType<Law> Basis, template<LawType, int, BasisType<Law>> typename T>
   struct law_of<T<Law,nDim,Basis>> : law_constant<Law> {};


// ---------- basis of variables ----------

/*
 * returns the spatial basis of a type, or false_type if no basis
 */
   template<typename T> struct basis_of; 

// remove cvrefs
   template<typename T> struct basis_of<const T> : basis_of<T> {};
   template<typename T> struct basis_of<volatile T> : basis_of<T> {};
   template<typename T> struct basis_of<const volatile T> : basis_of<T> {};

   template<typename T> struct basis_of<T&> : basis_of<T> {};
   template<typename T> struct basis_of<T&&> : basis_of<T> {};


   // helpers
   template<typename T> constexpr BasisType<law_of_v<T>> basis_of_v = basis_of<T>::value;
   template<typename T> using basis_of_t = typename basis_of<T>::type;
/*
 * returns whether a type has a basis or not
 */
   template<typename T>
   struct has_basis : std::conditional<std::is_same_v<typename basis_of<T>::type,
                                                      std::false_type>,
                                       std::false_type,
                                       std::true_type> {};

   // helper
   template<typename T> constexpr bool has_basis_v = has_basis<T>::value;

/*
 * returns whether multiple types have the same spatial dimension
 */
   template<typename T0, typename... Ts>
   using has_same_basis = all<std::is_same<basis_of_t<T0>,
                                           basis_of_t<Ts>>...>;

   // helper
   template<typename T0, typename... Ts>
   constexpr bool has_same_basis_v = has_same_basis<T0,Ts...>::value;

/*
 * specialisations of basis_of to extract basis from variable types
 */

   template<typename T>
   struct basis_of : std::false_type {};

   template<typename BasisType, BasisType Basis>
   struct basis_of<std::integral_constant<BasisType,Basis>> : std::integral_constant<BasisType,Basis> {};

   template<LawType Law, int nDim, BasisType<Law> Basis, template<LawType, int, BasisType<Law>> typename T>
   struct basis_of<T<Law,nDim,Basis>> : std::integral_constant<BasisType<Law>,Basis> {};


// ---------- return types of same law and spatial dimension ----------

   template<typename T> using    species_t =    Species< law_of_v<T> >;
   template<typename T> using      state_t =      State< law_of_v<T>, dim_of_v<T> >;
   template<typename T> using fluxresult_t = FluxResult< law_of_v<T>, dim_of_v<T> >;


// ---------- is type an instantiation of a particular type template? ----------

// Species
   template<typename T>
   struct is_Species : std::false_type {};

   template<LawType Law>
   struct is_Species<Species<Law>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_Species_v = is_Species<T>::value;

// State
   template<typename T>
   struct is_State : std::false_type {};

   template<LawType Law, int nDim>
   struct is_State<State<Law,nDim>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_State_v = is_State<T>::value;

// FluxResult
   template<typename T>
   struct is_FluxResult : std::false_type {};

   template<LawType Law, int nDim>
   struct is_FluxResult<FluxResult<Law,nDim>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_FluxResult_v = is_FluxResult<T>::value;

// VariableSet
   template<typename T>
   struct is_VariableSet : std::false_type {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct is_VariableSet<VariableSet<Law,nDim,Basis>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_VariableSet_v = is_VariableSet<T>::value;

// VariableDelta
   template<typename T>
   struct is_VariableDelta : std::false_type {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct is_VariableDelta<VariableDelta<Law,nDim,Basis>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_VariableDelta_v = is_VariableDelta<T>::value;

// VariableSet with specific Law and Basis
   template<typename T, LawType Law, BasisType<Law> Basis>
   struct is_specialised_VarSet : std::false_type {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct is_specialised_VarSet<VariableSet<Law,nDim,Basis>,Law,Basis> : std::true_type{};

   // helper
   template<typename T, LawType Law, BasisType<Law> Basis>
   constexpr bool is_specialised_VarSet_v = is_specialised_VarSet<T,Law,Basis>::value;

// VariableDelta with specific Law and Basis
   template<typename T, LawType Law, BasisType<Law> Basis>
   struct is_specialised_VarDelta : std::false_type {};

   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct is_specialised_VarDelta<VariableDelta<Law,nDim,Basis>,Law,Basis> : std::true_type{};

   // helper
   template<typename T, LawType Law, BasisType<Law> Basis>
   constexpr bool is_specialised_VarDelta_v = is_specialised_VarDelta<T,Law,Basis>::value;

