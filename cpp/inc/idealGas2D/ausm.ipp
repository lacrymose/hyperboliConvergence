# include <math.h>

namespace IdealGas2D
{
      inline void ausm( Species& gas, float n[3], ViscousVariables& ql, ViscousVariables& qr, ConservedVariables& f, float lmax )
  {
   // left/right states
      double ul,vl;
      double ur,vr;
      double unl, rl,pl,tl, al,hl,kl;
      double unr, rr,pr,tr, ar,hr,kr;

   // left/right numerical states/splittings
      double als,ml,mlp,plp;
      double ars,mr,mrm,prm;
      double m2m,m2p;

   // interface values
      double as,m2,m0,fa, ra,ma,pa;
      double psi[5];
      double pu,mp, delu,delp;

   // parameters/constants
      double alpha0,alpha,beta, Ku,Kp, sigma;
      double ascoeff;
      float  gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      alpha0=0.1875;
      beta  =0.125;
      Ku    =0.75;
      Kp    =0.25;
      sigma =1.0;

      ascoeff = 2.*(gas.gamma-1.)/(gas.gamma+1.);

   // left state

      // thermodynamic quantities
      ul = ql[0];
      vl = ql[1];
      tl = ql[2];
      pl = ql[3];

      rl = gas.Rgas*tl;
      al = rl*gas.gamma;
      rl = pl/rl;
      hl = al*gam1;
      al = sqrt(al);

      // kinetic energy
      kl = ul*ul;
      kl+= vl*vl;
      kl*= 0.5;
      hl+= kl;

      // velocities
      unl = n[0]*ul;
      unl+= n[1]*vl;

      als = ascoeff*hl;
      als = als / fmax( unl,sqrt(als) );

   // right state

      // thermodynamic quantities
      ur = qr[0];
      vr = qr[1];
      tr = qr[2];
      pr = qr[3];

      rr = gas.Rgas*tr;
      ar = rr*gas.gamma;
      rr = pr/rr;
      hr = ar*gam1;
      ar = sqrt(ar);

   // kinetic energy
      kr = ur*ur;
      kr+= vr*vr;
      kr*= 0.5;
      hr+= kr;

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
      fa = 1.;
      if( m0<1. ){ fa = m0*( 2. - m0 ); }

      alpha = alpha0*( 5.*fa*fa - 4. );

   // fourth order mach splitting, 5th order pressure splitting

   // left mach splitting
      if( ml < -1. )
     {
         mlp = 0.;
         plp = 0.;
         m2p = 0.;
         mrm = 0.;
     }
      else
     {
         if( ml < 1. )
        {
            m2p = ml+1.;
            m2p = 0.25*m2p*m2p;

            m2m = ml-1.;
            m2m =-0.25*m2m*m2m;

            mlp = m2p*( 1. - 16.*beta*m2m );
            plp = m2p*( 2. - ml - 16.*alpha*m2m*ml );
        }
         else
        {
            mlp = ml;
            plp = 1.;
            m2p = 0.;
            mrm = 0.;
        }
     }

   // right mach splitting
      if( mr < -1. )
     {
         mrm = mr;
         prm = 1.;
         m2p = 0.;
         mrm = 0.;
     }
      else
     {
         if( mr< 1. )
        {
            m2p = mr+1;
            m2p = 0.25*m2p*m2p;

            m2m = mr-1;
            m2m =-0.25*m2m*m2m;

            mrm = m2m*( 1. + 16.*beta*m2p );
            prm = m2m*(-2. - mr + 16.*alpha*m2p*mr );
        }
         else
        {
            mrm = 0.;
            prm = 0.;
            m2p = 0.;
            mrm = 0.;
        }
     }

   // mach number diffusion terms
      ra = 0.5*( rl + rr );
      delu = unr-unl;
      delp = pr - pl;

      mp = -Kp*fmax( 1.-sigma*m2, 0. );
      mp*=  delp/( fa*ra*as*as );

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

      double mdot = as*psi[0]*ma;
      mdot*=n[2];
      pa  *=n[2];

      f[0] = mdot;
      f[1] = psi[1]*mdot + n[0]*pa;
      f[2] = psi[2]*mdot + n[1]*pa;
      f[3] = psi[3]*mdot;

      lmax = fmax(         as,  fmax(       al,        ar  ) );
      lmax+= fmax( fabs(ma*as), fmax( fabs(unl), fabs(unr) ) );
      lmax*= n[2];

      return;
  }
}
