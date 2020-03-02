
namespace IdealGas2D
{
   inline void Slau::operator()( const Species& gas, const Types::Real n[3], const State& sl, const State& sr, ConservedDelta& f, Types::Real& lmax ) const
  {
      ConservedDelta  psil,psir;
      ConservedDelta      fc,fp;

      Types::Real     rl,ul,vl,pl,hl,al,unl,ml;
      Types::Real     rr,ur,vr,pr,hr,ar,unr,mr;

      Types::Real                           dr,dp;
      Types::Real                p_c,p_dp,p_du,pa;
      Types::Real               betal,betar,chi,g;
      Types::Real         una,ua,aa,aa1,ma,ma_hat;
      Types::Real     mdot_c,mdot_dr,mdot_dp,mdot;

      int  mlsup,mlsub;
      int  mrsup,mrsub;

   // left state
      rl = sl.density();
      ul = sl.velocityX();
      vl = sl.velocityY();
      pl = sl.pressure();
      hl = sl.specificTotalEnthalpy();
      al = sl.speedOfSound2();
      al = sqrt(al);

   // left state
      rr = sr.density();
      ur = sr.velocityX();
      vr = sr.velocityY();
      pr = sr.pressure();
      hr = sr.specificTotalEnthalpy();
      ar = sr.speedOfSound2();
      ar = sqrt(ar);

   // deltas
      dr = rr-rl;
      dp = pr-pl;

   // average velocities
      ua = 0.5*( ul*ul + vl*vl + ur*ur + vr*vr );
      ua = sqrt( ua );

      aa = 0.5*( al + ar );
      aa1= 1.0/aa;

      ma = ua*aa1;
      ma_hat = fmin( 1.0,ma );

   // left/right mach numbers
      unl = ul*n[0] + vl*n[1];
      unr = ur*n[0] + vr*n[1];

      ml = unl*aa1;
      mr = unr*aa1;

      una = ( rl*fabs(unl) + rr*fabs(unr) )/( rl+rr );

   // scaling parameters
      chi= ( 1.-ma_hat )*( 1.-ma_hat );

      g = -fmax(fmin(ml,0.),-1.)*fmin(fmax(mr,0.),1.);

   // interface mass flux
      mdot_c  = rl*unl + rr*unr;
      mdot_dr = una*dr;
      mdot_dp = chi*dp*aa1;

      mdot = 0.5*( ( mdot_c - mdot_dr )*(1.-g) - mdot_dp );

   // pressure splitting
      mlsub = std::signbit( ml-1 );
      mlsup = 1-mlsub;

      mrsub = std::signbit( mr-1 );
      mrsup = 1-mrsub;

      betal = mlsub*MachSplitting::P3(  1, ml );
      betal+= mlsup*MachSplitting::P1(  1, ml );

      betar = mrsub*MachSplitting::P3( -1, mr );
      betar+= mrsup*MachSplitting::P1( -1, mr );

   // interface pressure
      p_c  = ( pl + pr );
      p_dp = ( betal - betar )*dp;
//    p_du =-(1. - chi )*( betal + betar - 1. )*( pl + pr );  // slau
      p_du =         -ma*( betal + betar - 1. )*( pl + pr );  // slau2

      pa = 0.5*( p_c - p_dp - p_du );

   // scale flux by face area
      mdot*=n[2];
      pa  *=n[2];

   // convected quantities
      psil[0] = 1.;
      psil[1] = ul;
      psil[2] = vl;
      psil[3] = hl;

      psir[0] = 1.;
      psir[1] = ur;
      psir[2] = vr;
      psir[3] = hr;

   // convective flux
      fc = MachSplitting::M1(  1,mdot )*psil
         + MachSplitting::M1( -1,mdot )*psir;

   // pressure flux
      fp[0] = 0.;
      fp[1] = pa*n[0];
      fp[2] = pa*n[1];
      fp[3] = 0.;

   // assemble flux
      f = fc + fp;

   // spectral radius
      lmax = fmax(  aa, fmax(       al,        ar  ) );
      lmax+= fmax( una, fmax( fabs(unl), fabs(unr) ) );
      lmax*= n[2];
      lmax*=  2.0;

      return;
  }
}
