
# include <math.h>

namespace IdealGas2D
{
   inline void laxFriedrichs( Species& gas, float n[3], State& sl, State& sr, ConservedVariables& f, float& lmax )
  {
      ConservedVariables   fl,fr,fd,fc;
      ConservedVariables   ql,qr;
      float  ll,lr;

   // central flux
      exactFlux( gas, n, sl, fl, ll );
      exactFlux( gas, n, sr, fr, lr );

      lmax = fmax( ll, lr );

      fc = 0.5*( fl + fr );

   // upwind flux
      ql = conservedVariables( gas, sl );
      qr = conservedVariables( gas, sr );

      fd = -0.5*lmax*( qr - ql );

   // construct flux
      f = fc + fd;

      f*= n[2];

      return;
  }
}
