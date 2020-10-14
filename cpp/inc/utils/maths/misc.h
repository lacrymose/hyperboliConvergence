
# pragma once

# include <utils/concepts.h>

# include <array>
# include <cmath>

namespace utils
{

   template<typename ElemT, size_t N0>
   using matrix_1 = std::array<ElemT,N0>;

   template<typename ElemT, size_t N0, size_t N1>
   using matrix_2 = std::array<
                    std::array<ElemT,N1>,N0>;

   template<typename ElemT, size_t N0, size_t N1, size_t N2>
   using matrix_3 = std::array<
                    std::array<
                    std::array<ElemT,N2>,N1>,N0>;

   template<typename ElemT, size_t N0, size_t N1, size_t N2, size_t N3>
   using matrix_4 = std::array<
                    std::array<
                    std::array<
                    std::array<ElemT,N3>,N2>,N1>,N0>;


   template<floating_point Real>
   Real gaussian( Real mean, Real stddev, Real x )
  {
      constexpr Real sqrt2pi1 = 1./std::sqrt( 2.*M_PI );
      const Real stddev1 = 1./stddev;
      const Real pwr =  -0.5*( x-mean )*( x-mean )*stddev1*stddev1;
      return stddev1*sqrt2pi1*std::exp( pwr );
  }
}
