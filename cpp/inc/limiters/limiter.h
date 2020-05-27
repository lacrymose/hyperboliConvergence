
# pragma once

# include <utils/concepts.h>

namespace Limiters
{
/*
 * Piecewise Constant Method - 1st order extrapolation
 */
   struct NoLimit1
  {
      template<floating_point Real>
      Real operator()( const Real a, const Real b ) const;
  };

/*
 * Central Difference - 2nd order extrapolation
 */
   struct NoLimit2
  {
      template<floating_point Real>
      Real operator()( const Real a, const Real b ) const;
  };

/*
 * Parabolic reconstruction - 3rd order extrapolation
 */
   struct NoLimit3
  {
      template<floating_point Real>
      Real operator()( const Real a, const Real b ) const;
  };

/*
 * Van Albada's 2nd order limiter
 */
   struct VanAlbada2
  {
      template<floating_point Real>
      Real operator()( const Real a, const Real b ) const;
  };

/*
 * MinMod 2nd order limiter
 */
   struct MinMod2
  {
      template<floating_point Real>
      Real operator()( const Real a, const Real b ) const;
  };

/*
 * Van Leer's 2nd order limiter
 */
   struct VanLeer2
  {
      template<floating_point Real>
      Real operator()( const Real a, const Real b ) const;
  };

/*
 * Cada's 3rd order limiter
 */
   struct Cada3
  {
      template<floating_point Real>
      Real operator()( const Real a, const Real b ) const;
  };
}

# include <limiters/limiter.ipp>

