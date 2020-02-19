
namespace IdealGas2D
{
   inline void laxFriedrichs( const Species& gas, const float n[3], const State& sl, const State& sr, VariableSet<'c'>& f, float& lmax )
  {
      VariableSet<'c'>   fl,fr,fd,fc;
      VariableSet<'c'>   ql,qr;
      float  ll,lr;

   // central flux
      exactFlux( gas, n, sl, fl, ll );
      exactFlux( gas, n, sr, fr, lr );

      lmax = fmax( ll, lr );

      fc = 0.5*( fl + fr );

   // upwind flux
      ql = VariableSet<'c'>( gas, sl );
      qr = VariableSet<'c'>( gas, sr );

      fd = -0.5*lmax*( qr - ql );

   // construct flux
      f = fc + fd;

      f*= n[2];

      return;
  }
}
