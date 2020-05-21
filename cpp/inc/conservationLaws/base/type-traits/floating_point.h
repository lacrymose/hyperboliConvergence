
# pragma once

# include <conservationLaws/base/declarations.h>

# include <utils/concepts.h>
# include <utils/type-traits.h>

# include <type_traits>


// ---------- floating point precision of types ----------

/*
 * returns the spatial law of a type, or -1 if it has no law parameter
 */
   template<typename T> struct fptype_of;

// remove cvrefs
   template<typename T> struct fptype_of<const T> : fptype_of<T> {};
   template<typename T> struct fptype_of<volatile T> : fptype_of<T> {};
   template<typename T> struct fptype_of<const volatile T> : fptype_of<T> {};

   template<typename T> struct fptype_of<T&> : fptype_of<T> {};
   template<typename T> struct fptype_of<T&&> : fptype_of<T> {};


   // helpers
   template<typename T> using fptype_of_t = typename fptype_of<T>::type;


/*
 * returns whether a type has a spatial law or not
 */
   template<typename T>
   struct has_fptype : std::conditional_t<std::is_same_v<fptype_of_t<T>,
                                                         std::false_type>,
                                          std::false_type,
                                          std::true_type> {};

   // helpers
   template<typename T> constexpr bool has_fptype_v = has_fptype<T>::value;
// template<typename T> using has_fptype_t = typename has_fptype<T>::type;
   

/*
 * returns whether multiple types use the same floating point precision
 */
   template<typename T0, typename... Ts>
   using has_same_fptype = all< std::is_same<fptype_of_t<T0>,
                                             fptype_of_t<Ts>>... >;

   // helpers
   template<typename T0, typename...Ts>
   constexpr bool has_same_fptype_v = has_same_fptype<T0,Ts...>::value;


/*
 * specialisations of fptype_of to extract fptype from template types with typical parameter lists
 */

   // fptype_of is false if template argument does not have a floating point precision
// template<typename T>
// struct fptype_of : std::false_type {};

   template<typename T>
   struct fptype_of : std::conditional_t< std::is_floating_point_v<T>,
                                          type_identity<T>,
                                          std::false_type> {};

   // covers all current use cases, and probably many future ones?
   template<LawType Law, floating_point Real, template<LawType, floating_point> typename T>
   struct fptype_of<T<Law,Real>> : type_identity<Real> {};

   template<LawType Law, int nDim, floating_point Real, template<LawType, int, floating_point> typename T>
   struct fptype_of<T<Law,nDim,Real>> : type_identity<Real> {};

   template<LawType Law, int nDim, BasisType<Law> Basis, floating_point Real, template<LawType, int, BasisType<Law>, floating_point> typename T>
   struct fptype_of<T<Law,nDim,Basis,Real>> : type_identity<Real> {};



