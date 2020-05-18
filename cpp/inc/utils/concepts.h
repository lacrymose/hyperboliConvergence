
# pragma once

# include <utils/type-traits.h>

/*
 * The concept same_template_as<T,U> is satisfied if and only if T and U denote the same type template
 *    shadows std::same_as<T,U> concept, but for type templates not types
 */
   template<template<typename...> typename T, template<typename...> typename U>
   concept bool same_template_as = is_same_template_v<T,U> && is_same_template_v<U,T>;

/*
 * checks all types in pack are the same
 */
// template<typename T, typename... Us>
// concept bool Same = all< std::is_same<T,Us>... >::value;

   template<typename T, typename... Us>
   concept bool Same = ( std::is_same_v<T,Us> && ... );

