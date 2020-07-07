
# pragma once

# include <conservationLaws/base/base.h>
# include <conservationLaws/base/concepts.h>

# include <utils/concepts.h>

namespace Limiters
{
/*
 * Base class for CRTP to give an interface to extrapolate VariableDeltas
 */
   template<typename Limiter>
   struct VariableLimiterInterface
  {
      template<ImplementedVarDelta VarDelT>
      VarDelT operator()( const VarDelT& dq0, const VarDelT& dq1 ) const
     {
         constexpr int N=VarDelT::N;
         VarDelT dqlim;
         for( int i=0; i<N; i++ )
        {
            dqlim[i] = static_cast<const Limiter*>(this)->limit( dq0[i],dq1[i] );
        }
         return dqlim;
     }
  };

   template<typename Limiter>
   struct ScalarLimiterInterface
  {
   // g++ says this call is ambiguous with the VariableLimiterInterface operator()
   // I don't know why yet. need to investigate.
//    template<floating_point Real>
//    float operator()( const float a, const float b ) const
//   {
//       return static_cast<const Limiter*>(this)->limit( a,b );
//   }
  };


/*
 * Piecewise Constant Method - 1st order extrapolation
 */
   struct NoLimit1 : ScalarLimiterInterface<  NoLimit1>,
                     VariableLimiterInterface<NoLimit1>
  {
      template<floating_point Real>
      Real limit( const Real a, const Real b ) const;
  };

/*
 * Central Difference - 2nd order extrapolation
 */
   struct NoLimit2 : ScalarLimiterInterface<  NoLimit2>,
                     VariableLimiterInterface<NoLimit2>
  {
      template<floating_point Real>
      Real limit( const Real a, const Real b ) const;
  };

/*
 * Parabolic reconstruction - 3rd order extrapolation
 */
   struct NoLimit3 : ScalarLimiterInterface<  NoLimit3>,
                     VariableLimiterInterface<NoLimit3>
  {
      template<floating_point Real>
      Real limit( const Real a, const Real b ) const;
  };

/*
 * Van Albada's 2nd order limiter
 */
   struct VanAlbada2 : ScalarLimiterInterface<  VanAlbada2>,
                       VariableLimiterInterface<VanAlbada2>
  {
      template<floating_point Real>
      Real limit( const Real a, const Real b ) const;
  };

/*
 * MinMod 2nd order limiter
 */
   struct MinMod2 : ScalarLimiterInterface<  MinMod2>,
                    VariableLimiterInterface<MinMod2>
  {
      template<floating_point Real>
      Real limit( const Real a, const Real b ) const;
  };

/*
 * Superbee 2nd order limiter
 */
   struct Superbee2 : ScalarLimiterInterface<  Superbee2>,
                      VariableLimiterInterface<Superbee2>
  {
      template<floating_point Real>
      Real limit( const Real a, const Real b ) const;
  };

/*
 * Monotonized Central 2nd order limiter
 */
   struct MonotonizedCentral2 : ScalarLimiterInterface<  MonotonizedCentral2>,
                                VariableLimiterInterface<MonotonizedCentral2>
  {
      template<floating_point Real>
      Real limit( const Real a, const Real b ) const;
  };

/*
 * Van Leer's 2nd order limiter
 */
   struct VanLeer2 : ScalarLimiterInterface<  VanLeer2>,
                     VariableLimiterInterface<VanLeer2>
  {
      template<floating_point Real>
      Real limit( const Real a, const Real b ) const;
  };

/*
 * Cada's 3rd order limiter
 */
   struct Cada3 : ScalarLimiterInterface<  Cada3>,
                  VariableLimiterInterface<Cada3>
  {
      template<floating_point Real>
      Real limit( const Real a, const Real b ) const;
  };
}

# include <limiters/limiter.ipp>

