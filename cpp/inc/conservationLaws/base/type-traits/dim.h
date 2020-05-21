
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


