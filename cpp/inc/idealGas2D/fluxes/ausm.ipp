
namespace IdealGas2D
{
   inline void Ausm::operator()( const Species& gas, const Types::Real n[3], const State& sl, const State& sr, ConservedDelta& f, Types::Real& lmax ) const
  {
   // left/right states
      Types::Real ul,vl;
      Types::Real ur,vr;
      Types::Real unl, rl,pl, al,hl;
      Types::Real unr, rr,pr, ar,hr;

   // left/right numerical states/splittings
      Types::Real als,ml,mlp,plp;
      Types::Real ars,mr,mrm,prm;

   // interface values
      Types::Real as,m2,m0,fa, ra,ma,pa;
      Types::Real psi[5];
      Types::Real pu,mp, delu,delp;
      Types::Real alpha;

   // parameters/constants
      Types::Real ascoeff;
      Types::Real  gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      ascoeff = 2.*(gas.gamma-1.)/(gas.gamma+1.);

   // left state

      ul = sl.velocityX();
      vl = sl.velocityY();
      pl = sl.pressure();
      rl = sl.density();
      hl = sl.specificTotalEnthalpy();
      al = sqrt( sl.speedOfSound2() );

      // velocities
      unl = n[0]*ul;
      unl+= n[1]*vl;

      als = ascoeff*hl;
      als = als / fmax( unl,sqrt(als) );

   // right state

      ur = sr.velocityX();
      vr = sr.velocityY();
      pr = sr.pressure();
      rr = sr.density();
      hr = sr.specificTotalEnthalpy();
      ar = sqrt( sr.speedOfSound2() );

      // velocities
      unr = n[0]*ur;
      unr+= n[1]*vr;

      ars = ascoeff*hr;
      ars = ars / fmax( -unr,sqrt(ars) );

   // interface mach number splitting
      as = fmin( ars, als );

      ml = unl/as;
      mr = unr/as;

      m2 = 0.5*( ml*ml + mr*mr );
      m0 = fmin( 1., fmax( m2, gas.minf*gas.minf ) );
      m0 = sqrt( m0 );
      fa = m0*( 2. - m0 );

      alpha = alpha0*( 5.*fa*fa - 4. );

   // fourth order mach splitting, 5th order pressure splitting
      int mlsup,mlsub;
      int mrsup,mrsub;

      mlsup = fabs(ml)>1 ? 1:0;
      mlsub = 1-mlsup;

      mrsup = fabs(mr)>1 ? 1:0;
      mrsub = 1-mrsup;

   // mach splittings
      mlp = mlsub*MachSplitting::M4(  1, ml, beta );
      mlp+= mlsup*MachSplitting::M1(  1, ml );

      mrm = mrsub*MachSplitting::M4( -1, mr, beta );
      mrm+= mrsup*MachSplitting::M1( -1, mr );

   // pressure splittings
      plp = mlsub*MachSplitting::P5(  1, ml, alpha );
      plp+= mlsup*MachSplitting::P1(  1, ml );

      prm = mrsub*MachSplitting::P5( -1, mr, alpha );
      prm+= mrsup*MachSplitting::P1( -1, mr );

   // mach number diffusion terms
      ra = 0.5*( rl + rr );
      delu = unr-unl;
      delp = pr - pl;

      mp = -Kp*fmax( 1.-sigma*m2, 0. );
      mp*=          delp/( fa*ra*as*as );
//    mp*=  (1.-fa)*delp/(    ra*as*as );   // cfl limit bodge -not tested for proper low-mach scaling!

      pu = -Ku*plp*prm;
      pu*=  2.*ra*(fa*as)*delu;

   // interface mach/pressure
      ma = mlp    + mrm    + mp;
      pa = plp*pl + prm*pr + pu;

      if( ma > 0 )
     {
         psi[0]=rl;
         psi[1]=ul;
         psi[2]=vl;
         psi[3]=hl;
     }
      else
     {
         psi[0] = rr;
         psi[1] = ur;
         psi[2] = vr;
         psi[3] = hr;
     }

      Types::Real mdot = as*psi[0]*ma;
      mdot*=n[2];
      pa  *=n[2];

      f[0] = mdot;
      f[1] = psi[1]*mdot + n[0]*pa;
      f[2] = psi[2]*mdot + n[1]*pa;
      f[3] = psi[3]*mdot;

      Types::Real  lmd=1.;
      lmd= 0.5*(m0+1)/fa;

      lmax = lmd*fmax(     as,  fmax(       al,        ar  ) );
      lmax+= fmax( fabs(ma*as), fmax( fabs(unl), fabs(unr) ) );
      lmax*= n[2];

      return;
  }
}
