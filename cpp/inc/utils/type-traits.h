
# pragma once

# include <type_traits>

# include <cstddef>

// type trait utilities

// -------- returns type of argument --------

   // to be replaced by std::type_identity in C++20
   template<typename T>
   struct type_identity { using type=T; };

   // helper
   template<typename T>
   using type_identity_t = typename type_identity<T>::type;

// -------- returns whether type is an integer representaion --------

   template<typename T>
   struct is_integer : std::false_type {};

   template<> struct is_integer<          int> : std::true_type {};
   template<> struct is_integer<    short int> : std::true_type {};
   template<> struct is_integer<     long int> : std::true_type {};
   template<> struct is_integer<long long int> : std::true_type {};

   template<> struct is_integer<unsigned           int> : std::true_type {};
   template<> struct is_integer<unsigned     short int> : std::true_type {};
   template<> struct is_integer<unsigned      long int> : std::true_type {};
   template<> struct is_integer<unsigned long long int> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_integer_v = is_integer<T>::value;

// -------- remove const/volative qualifiers and references --------

   // these will be replaced by std::remove_cvref in C++20
   template<typename T>
   struct remove_cvref : std::remove_cv<std::remove_reference_t<T>> {};

   template<typename T> using remove_cvref_t = typename remove_cvref<T>::type;

// -------- test if two type templates are the same --------

   // mirrors std::is_same<T,U> but for type templates not types
   template<template<typename...> typename T, template<typename...> typename U>
   struct is_same_template : std::false_type {};

   template<template<typename...> typename T>
   struct is_same_template<T,T> : std::true_type {};

   // helper
   template<template<typename...> typename T, template<typename...> typename U>
   constexpr bool is_same_template_v = is_same_template<T,U>::value;

// -------- test if value member is of bool type --------

   template<typename T>
   constexpr bool has_bool_value_v = std::is_same_v<typename T::value_type,bool>;


// -------- test pack size --------

   template<std::size_t N, typename... Ts>
   struct is_pack_of_n : std::bool_constant<sizeof...(Ts)==N> {};

   template<typename... Ts> using is_empty_pack  = is_pack_of_n<0,Ts...>;
   template<typename... Ts> using is_pack_of_one = is_pack_of_n<1,Ts...>;


// -------- all of conditions are true --------

   template<typename... Conds>
   struct all : std::true_type {};

   template<typename Cond0, typename... Conds>
      requires has_bool_value_v<Cond0>
   struct all<Cond0,Conds...> : std::conditional_t< Cond0::value,
                                                    all<Conds...>,
                                                    std::false_type > {};


// -------- none of conditions are true --------

   template<typename... Conds>
   struct none : std::true_type {};

   template<typename Cond0, typename... Conds>
      requires has_bool_value_v<Cond0>
   struct none<Cond0,Conds...> : std::conditional_t< Cond0::value,
                                                     std::false_type,
                                                     none<Conds...> > {};


// -------- any of conditions are true --------

   template<typename... Conds>
   struct any : std::false_type {};

   template<typename Cond0, typename... Conds>
      requires has_bool_value_v<Cond0>
   struct any<Cond0,Conds...> : std::conditional_t< Cond0::value,
                                                    std::true_type,
                                                    any<Conds...> > {};

