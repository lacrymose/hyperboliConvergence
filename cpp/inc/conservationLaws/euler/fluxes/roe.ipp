
   template<int nDim, floating_point Real>
   State<LawType::Euler,nDim,Real> roeAverage( const Species<LawType::Euler,Real>& species,
                                               const State<LawType::Euler,nDim,Real>&   sl,
                                               const State<LawType::Euler,nDim,Real>&   sr )
  {
      State<LawType::Euler,nDim,Real> ravg;


   // left/right weighting functions
      const Real cr = sqrt( sr.density()/sl.density() );
      const Real wl = 1./( 1.+cr );
      const Real wr = 1.-wl;

   // roe average density
      const Real ra = cr*sl.density();

   // roe average velocities
      Real ka=0;
      for( unsigned int i=0; i<nDim; i++ )
     {
         const Real ua = wl*sl.velocity(i)
                        +wr*sr.velocity(i);
         ka+=ua*ua;
         ravg.velocity(i)=ua;
     }
      ka*=0.5;

   // roe average enthalpy
      const Real ha = wl*sl.specificTotalEnthalpy()
                     +wr*sr.specificTotalEnthalpy();

   // roe average speed of sound(^2)
      const Real a2a = (species.gamma-1.)*( ha-ka );

   // roe average pressure and temperature
      const Real ta = a2a/(species.R*species.gamma);
      const Real pa = ra*species.R*ta;

   // assemble
      ravg.pressure()    = pa;
      ravg.density()     = ra;
      ravg.temperature() = ta;
      ravg.specificTotalEnthalpy() = ha;
      ravg.velocity2()   = 2.*ka;
      ravg.speedOfSound2() = a2a;

      return ravg;
  }

   template<int nDim, floating_point Real>
   WaveSpeeds<LawType::Euler,nDim,Real> entropyfix( const WaveSpeeds<LawType::Euler,nDim,Real>& la,
                                                    const WaveSpeeds<LawType::Euler,nDim,Real>& ll,
                                                    const WaveSpeeds<LawType::Euler,nDim,Real>& lr )
  {
      WaveSpeeds<LawType::Euler,nDim,Real> speeds;

   // no fix needed for convective speeds
      const Real ua = fabs(la[0]);
      speeds[0]=ua;
      for( unsigned int i=0; i<nDim-1; i++ ){ speeds[3+i] = ua; }

   // harten's fix for nonlinear speeds
      speeds[1]=fabs(la[1]);
      speeds[2]=fabs(la[2]);

      constexpr Real eps=0.05;
      for( unsigned int i=1; i<3; i++ )
     {
         const Real le = eps*fmax(  0.,
                              fmax( la[i]-ll[i],
                                    lr[i]-la[i] ) );

         const Real laa = fabs(la[i]);
         speeds[i] = (laa>le) ? laa : 0.5*( le + laa*laa/le );
     }

      return speeds;
  }

