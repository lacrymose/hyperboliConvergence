
# pragma once

# include <conservationLaws/base/declarations.h>

# include <utils/type-traits.h>

# include <type_traits>


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



