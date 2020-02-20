
namespace IdealGas2D
{
   inline void LaxFriedrichs::operator()( const Species& gas, const float n[3], const State& sl, const State& sr, ConservedDelta& f, float& lmax ) const
  {
      ConservedDelta   fl,fr,fd,fc;
      ConservedVariables     ql,qr;
      float  ll,lr;

   // central flux
      exactFlux( gas, n, sl, fl, ll );
      exactFlux( gas, n, sr, fr, lr );

      lmax = fmax( ll, lr );

      fc = 0.5*( fl + fr );

   // upwind flux
      ql = ConservedVariables( gas, sl );
      qr = ConservedVariables( gas, sr );

      fd = -0.5*lmax*( qr - ql );

   // construct flux
      f = fc + fd;

      f*= n[2];

      return;
  }
}
