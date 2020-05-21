
# pragma once

# include <conservationLaws/base/declarations.h>

# include <utils/type-traits.h>

# include <type_traits>


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

   template<LawType Law, int nDim, BasisType<Law> Basis, floating_point Real, template<LawType, int, BasisType<Law>, floating_point> typename T>
   struct basis_of<T<Law,nDim,Basis,Real>> : std::integral_constant<BasisType<Law>,Basis> {};



