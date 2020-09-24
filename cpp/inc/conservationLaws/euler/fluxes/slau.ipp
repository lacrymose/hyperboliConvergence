
# include <conservationLaws/euler/fluxes/machSplittings.ipp>

      template<LowMachScaling VFluxScaling,
               LowMachScaling PFluxScaling>                       // struct template parameters

      template<EulerState StateT, int nDim, floating_point Real>  // member function template parameters
         requires   SameDim<   StateT,dim_constant<nDim>>
                 && SameFPType<StateT,Real>

      FluxResult<LawType::Euler,nDim,Real>    // return type
               Slau<VFluxScaling,PFluxScaling>::flux( const Species<LawType::Euler,Real>& species,
                                                      const geom::Surface<nDim,Real>&        face,
                                                      const StateT&                            sl,
                                                      const StateT&                            sr ) const
  {
// average velocities
      const Real ua = sqrt( 0.5*( sl.velocity2() + sr.velocity2() ) );

      const Real aa = 0.5*(  sqrt( sl.speedOfSound2() )
                           + sqrt( sr.speedOfSound2() ) );
      const Real aa1= 1.0/aa;

      const Real ma = ua*aa1;
      const Real ma_hat = fmin( 1.0,ma );

// left/right mach numbers
      const Real unl = projectedVelocity( face, sl );
      const Real unr = projectedVelocity( face, sr );

      const Real ml = unl*aa1;
      const Real mr = unr*aa1;

      const Real una = ( sl.density()*fabs(unl) + sr.density()*fabs(unr) )
                      /( sl.density()           + sr.density()           );

// scaling parameters
      const Real chi= ( 1.-ma_hat )*( 1.-ma_hat );

   // cut-off mach number
      const Real minf = species.minf;

   // ratio of acoustic timescale to simulation timestep
      const Real mu = species.lref/(aa*species.dt);

   // convective / adaptive  scaling parameters
      [[maybe_unused]] const Real m0_c = std::clamp<Real>(       ma,       minf, 1. );
      [[maybe_unused]] const Real m0_a = std::clamp<Real>( fmax( ma, mu ), minf, 1. );

   // scaling parameter for pressure diffusion in mass flux
      [[maybe_unused]] const Real m0_p = [&]() -> Real
     {
         if      constexpr( PFluxScaling == LowMachScaling::Convective ){ return m0_c; }
         else if constexpr( PFluxScaling == LowMachScaling::Acoustic   ){ return 1.0;  }
         else            /* PFluxScaling == LowMachScaling::Adaptive  */{ return m0_a; }
     }();

   // scaling parameter for velocity diffusion in pressure flux
      [[maybe_unused]] const Real m0_u = [&]() -> Real
     {
         if      constexpr( VFluxScaling == LowMachScaling::Convective ){ return m0_c; }
         else if constexpr( VFluxScaling == LowMachScaling::Acoustic   ){ return 1.0;  }
         else            /* VFluxScaling == LowMachScaling::Adaptive  */{ return m0_a; }
     }();

   // switch function to prevent negative density in very high speed expansions
      const Real g = -fmax(fmin(ml,0.),-1.)*fmin(fmax(mr,0.),1.);

// interface mass flux
      const Real dr = sr.density() - sl.density();
      const Real dp = sr.pressure()- sl.pressure();

   // central contribution
      const Real mdot_c  =  sl.density()*unl
                          + sr.density()*unr;

   // diffusive contributions
      const Real mdot_dr = una*dr;
   // const Real mdot_dr = 0.;

      const Real mdot_dp = chi*dp*aa1;       // slau
   // const Real mdot_dp = chi*dp/(aa*m0_p); // slau-p'
   // const Real mdot_dp = 0.;

   // assemble
      const Real mdot = 0.5*( ( mdot_c - mdot_dr )*(1.-g) - mdot_dp )*face.area;
   // const Real mdot = 0.5*mdot_c*face.area;

// interface pressure

   // pressure splitting
      const unsigned int mlsup = fabs(ml)>1 ? 1:0;
      const unsigned int mlsub = 1-mlsup;

      const unsigned int mrsup = fabs(mr)>1 ? 1:0;
      const unsigned int mrsub = 1-mrsup;

      const Real betal =  mlsub*MachSplit_P3(  1, ml )
                        + mlsup*MachSplit_P1(  1, ml );

      const Real betar =  mrsub*MachSplit_P3( -1, mr )
                        + mrsup*MachSplit_P1( -1, mr );

   // central component
      const Real p_c  = ( sl.pressure() + sr.pressure() );

   // diffusive contributions
      const Real p_dp = ( betal - betar )*dp;
   // const Real p_dp = 0.;

   // const Real p_du =-(1. - chi )*( betal + betar - 1. )*( p_c );  // slau
      const Real p_du =         -ma*( betal + betar - 1. )*( p_c );  // slau2
   // const Real p_du =       -m0_u*( betal + betar - 1. )*( p_c );  // slau-u'
   // const Real p_du = 0.;

   // assemble
      const Real pa = 0.5*( p_c - p_dp - p_du )*face.area;
   // const Real pa = 0.5*p_c*face.area;

// interface flux
      FluxResult<LawType::Euler,nDim,Real> result;

   // switch functions
      const int m1l = (ua>0) ? 1 : 0;
      const int m1r = (ua>0) ? 0 : 1;

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

   // zha-bilgen energy flux
//    const Real el = sl.specificTotalEnthalpy()-0.5*sl.velocity2();
//    const Real er = sr.specificTotalEnthalpy()-0.5*sr.velocity2();

//    result.flux[nDim+1] =( 0.5*(mdot_c-mdot_dr)*(  m1l*el
//                                                 + m1r*er )
//                         - 0.5*mdot_dp*(  m1l*sl.specificTotalEnthalpy()
//                                        + m1r*sr.specificTotalEnthalpy() )
//                         - 0.5*una*( p_dp + p_du ) )*face.area;

   // spectral radius
      const Real lmax =  sqrt( fmax( sl.speedOfSound2(),
                                     sr.speedOfSound2() ) )
                       + fmax( fabs(unl),
                               fabs(unr) );

      result.lambda = lmax*face.area;

      return result;
  }
