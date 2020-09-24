
# pragma once

# include <utils/concepts.h>

# include <array>
# include <cmath>

namespace utils
{

   template<typename ElemT, auto... Ns>
   struct matrixHelper;

   template<typename ElemT, size_t N0>
   struct matrixHelper<ElemT,N0>
  {
      using type = std::array<ElemT,N0>;
  };

   template<typename ElemT, size_t N0, size_t N1>
   struct matrixHelper<ElemT,N0,N1>
  {
      using type = std::array<
                   std::array<ElemT,N1>,N0>;
  };

   template<typename ElemT, size_t N0, size_t N1, size_t N2>
   struct matrixHelper<ElemT,N0,N1,N2>
  {
      using type = std::array<
                   std::array<
                   std::array<ElemT,N2>,N1>,N0>;
  };

   template<typename ElemT, size_t N0, size_t N1, size_t N2, size_t N3>
   struct matrixHelper<ElemT,N0,N1,N2,N3>
  {
      using type = std::array<
                   std::array<
                   std::array<
                   std::array<ElemT,N3>,N2>,N1>,N0>;
  };

   template<typename ElemT, auto... Ns>
   using matrix_t = typename matrixHelper<ElemT,Ns...>::type;


   template<floating_point Real>
   Real gaussian( Real mean, Real stddev, Real x )
  {
      constexpr Real sqrt2pi1 = 1./std::sqrt( 2.*M_PI );
      const Real stddev1 = 1./stddev;
      const Real pwr =  -0.5*( x-mean )*( x-mean )*stddev1*stddev1;
      return stddev1*sqrt2pi1*std::exp( pwr );
  }
}
