
# include <conservationLaws/euler/fluxes/machSplittings.ipp>

# include <algorithm>

      template<LowMachScaling VFluxScaling,                       // struct template parameters
               LowMachScaling PFluxScaling>

      template<EulerState StateT, int nDim, floating_point Real>  // member function template parameters
         requires   SameDim<   StateT,dim_constant<nDim>>
                 && SameFPType<StateT,Real>

      FluxResult<LawType::Euler,nDim,Real>    // return type
               AusmPlusUP<VFluxScaling,PFluxScaling>::flux( const Species<LawType::Euler,Real>& species,
                                                            const geom::Surface<nDim,Real>&        face,
                                                            const StateT&                            sl,
                                                            const StateT&                            sr ) const
  {
   // parameters/constants

      const Real alpha0 = 0.1875;
      const Real beta   = 0.125;
      const Real Ku     = 0.75*2.;
      const Real Kp     = 0.25;
      const Real sigma  = 1.0;

      const Real ascoeff = 2.*(species.gamma-1.)/(species.gamma+1.);

   // velocities

      const Real unl = projectedVelocity( face.metric[0], sl );
      const Real unr = projectedVelocity( face.metric[0], sr );

   // critical sonic speeds adjusted to prevent sonic expansion shock

      const Real als2 =  ascoeff*sl.specificTotalEnthalpy();
      const Real ars2 =  ascoeff*sr.specificTotalEnthalpy();

      const Real als  =  als2 / fmax(  unl,sqrt(als2) );
      const Real ars  =  ars2 / fmax( -unr,sqrt(ars2) );

   // interface speed of sound and density

      const Real ra = 0.5*(  sl.density()
                           + sr.density() );

      const Real as = fmin( ars, als );

   // mach numbers

      // left/right
      const Real ml = unl/as;
      const Real mr = unr/as;

      // cut-off mach number
      const Real minf2 = species.minf*species.minf;

      // local mach number
      const Real m2 = 0.5*( ml*ml + mr*mr );

      // ratio of acoustic timescale to simulation timestep
      const Real mu = species.lref/(as*species.dt);

      // convective / adaptive  scaling parameters
      [[maybe_unused]] const Real m0_c = sqrt( std::clamp<Real>(       m2,          minf2, 1. ) );
      [[maybe_unused]] const Real m0_a = sqrt( std::clamp<Real>( fmax( m2, mu*mu ), minf2, 1. ) );

   // scaling parameter for pressure diffusion in mass flux
      const Real m0_p = [&]() -> Real
     {
         if      constexpr( PFluxScaling == LowMachScaling::Convective ){ return m0_c; }
         else if constexpr( PFluxScaling == LowMachScaling::Acoustic   ){ return 1.0;  }
         else            /* PFluxScaling == LowMachScaling::Adaptive  */{ return m0_a; }
     }();

   // scaling parameter for velocity diffusion in pressure flux
      const Real m0_u = [&]() -> Real
     {
         if      constexpr( VFluxScaling == LowMachScaling::Convective ){ return m0_c; }
         else if constexpr( VFluxScaling == LowMachScaling::Acoustic   ){ return 1.0;  }
         else            /* VFluxScaling == LowMachScaling::Adaptive  */{ return m0_a; }
     }();

      const Real fa_p = m0_p*( 2. - m0_p );
      const Real fa_u = m0_u*( 2. - m0_u );

      const unsigned int mlsup = fabs(ml)>1 ? 1:0;
      const unsigned int mlsub = 1-mlsup;

      const unsigned int mrsup = fabs(mr)>1 ? 1:0;
      const unsigned int mrsub = 1-mrsup;

   // interface mach number

      // mach splittings (4th order)
      const Real mlp =  mlsub*MachSplit_M4(  1, ml, beta )
                      + mlsup*MachSplit_M1(  1, ml );

      const Real mrm =  mrsub*MachSplit_M4( -1, mr, beta )
                      + mrsup*MachSplit_M1( -1, mr );

      // diffusion terms
      const Real delp =  sr.pressure()
                       - sl.pressure();

      const Real mp = Kp * fmax( 1.-sigma*m2, 0. ) * delp / ( fa_p*ra*as*as );

      // assemble
      const Real ua = as*( mlp + mrm - mp )*face.area;

   // interface pressure

      // pressure splittings (5th order)
      const Real alpha = alpha0*( 5.*fa_u*fa_u - 4. );

      const Real plp =  mlsub*MachSplit_P5(  1, ml, alpha )
                      + mlsup*MachSplit_P1(  1, ml );

      const Real prm =  mrsub*MachSplit_P5( -1, mr, alpha )
                      + mrsup*MachSplit_P1( -1, mr );

      // diffusion terms
      const Real delu = unr-unl;
      const Real pu = Ku * plp*prm * ra*(fa_u*as) * delu;

      // assemble
      const Real pa = ( plp*sl.pressure() + prm*sr.pressure() - pu )*face.area;

   // interface flux
      FluxResult<LawType::Euler,nDim,Real> result;

   // switch functions
      const int m1l = (ua>0) ? 1 : 0;
      const int m1r = (ua>0) ? 0 : 1;

   // mass flux
      const Real mdot = ua*(  m1l*sl.density()
                            + m1r*sr.density() );

   // momentum fluxes
      for( int i=0; i<nDim; i++ )
     {
         result.flux[i] =  mdot*( m1l*sl.velocity(i)      // upwind velocity
                                 +m1r*sr.velocity(i) )
                         + pa*face.metric[0][i];          // pressure flux
     }

   // density flux
      result.flux[nDim] = mdot;

   // energy flux
      result.flux[nDim+1] = mdot*(  m1l*sl.specificTotalEnthalpy()
                                  + m1r*sr.specificTotalEnthalpy() );

//    spectral radius scaling for low mach numbers
      const Real lmd=1.;
//    const Real lmd= 0.5*(m0_p+1)/fa_p;

      const Real amax = sqrt( fmax(  as*as,
                               fmax( sl.speedOfSound2(),
                                     sr.speedOfSound2() ) ) );

      const Real lmax =  lmd*amax
                       + fmax(  fabs(ua),
                          fmax( fabs(unl),
                                fabs(unr) ) );

      result.lambda = lmax*face.area;

      return result;
  }
