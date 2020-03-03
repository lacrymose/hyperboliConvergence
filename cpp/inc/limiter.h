# ifndef LIMITER_H
# define LIMITER_H

# include <idealGas2D/idealGas2D.h>

# include <types.h>

# include <cmath>

namespace Limiters
{
/*
 * Limits left and right vector extrapolation across face
 */
   template<typename Limiter>
   struct FaceVectorLimiter
  {
      Limiter  limiter;
      template<typename DeltaType>
      inline void operator()( const IdealGas2D::VariableDelta<DeltaType>& dl,
                              const IdealGas2D::VariableDelta<DeltaType>& dc,
                              const IdealGas2D::VariableDelta<DeltaType>& dr,
                                    IdealGas2D::VariableDelta<DeltaType>& dl_lim,
                                    IdealGas2D::VariableDelta<DeltaType>& dr_lim ) const;
  };

/*
 * Limits left and right extrapolation across face
 */
   template<typename Limiter>
   struct FaceLimiter
  {
      Limiter  limiter;
      inline void operator()( Types::Real  dl,    Types::Real dc, Types::Real  dr,
                              Types::Real& dl_lim,                Types::Real& dr_lim ) const;
  };

/*
 * Piecewise Constant Method - 1st order extrapolation
 */
   struct NoLimit1
  {
      inline Types::Real operator()( Types::Real a, Types::Real b ) const;
  };

/*
 * Central Difference - 2nd order extrapolation
 */
   struct NoLimit2
  {
      inline Types::Real operator()( Types::Real a, Types::Real b ) const;
  };

/*
 * Parabolic reconstruction - 3rd order extrapolation
 */
   struct NoLimit3
  {
      inline Types::Real operator()( Types::Real a, Types::Real b ) const;
  };

/*
 * Van Albada's 2nd order limiter
 */
   struct VanAlbada2
  {
      inline Types::Real operator()( Types::Real a, Types::Real b ) const;
  };

/*
 * MinMod 2nd order limiter
 */
   struct MinMod2
  {
      inline Types::Real operator()( Types::Real a, Types::Real b ) const;
  };

/*
 * Van Leer's 2nd order limiter
 */
   struct VanLeer2
  {
      inline Types::Real operator()( Types::Real a, Types::Real b ) const;
  };

/*
 * Cada's 3rd order limiter
 */
   struct Cada3
  {
      Types::Real alpha;
      Types::Real  beta;

      inline Types::Real operator()( Types::Real a, Types::Real b ) const;
  };

}

# include <limiter.ipp>

# endif
