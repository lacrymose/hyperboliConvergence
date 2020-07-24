
# pragma once

# include<geometry/type-traits/dim.h>
# include<geometry/type-traits/floating_point.h>

# include <geometry/declarations.h>

# include <utils/type-traits.h>

# include <type_traits>

namespace geom
{
// ---------- is type an instantiation of a particular type template? ----------

// point

   template<typename T>
   struct is_Point : std::false_type {};

   template<int nDim, floating_point Real>
   struct is_Point<Point<nDim,Real>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_Point_v = is_Point<T>::value;

// direction

   template<typename T>
   struct is_Direction : std::false_type {};

   template<int nDim, floating_point Real>
   struct is_Direction<Direction<nDim,Real>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_Direction_v = is_Direction<T>::value;

// metric

   template<typename T>
   struct is_Metric : std::false_type {};

   template<int nDim, floating_point Real>
   struct is_Metric<Metric<nDim,Real>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_Metric_v = is_Metric<T>::value;

// surface

   template<typename T>
   struct is_Surface : std::false_type {};

   template<int nDim, floating_point Real>
   struct is_Surface<Surface<nDim,Real>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_Surface_v = is_Surface<T>::value;

// volume

   template<typename T>
   struct is_Volume : std::false_type {};

   template<int nDim, floating_point Real>
   struct is_Volume<Volume<nDim,Real>> : std::true_type {};

   // helper
   template<typename T>
   constexpr bool is_Volume_v = is_Volume<T>::value;

}
