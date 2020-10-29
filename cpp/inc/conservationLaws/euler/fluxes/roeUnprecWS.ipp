
   template<EulerState StateT, int nDim, floating_point Real>
      requires   SameDim<   StateT,dim_constant<nDim>>
              && SameFPType<StateT,Real>
   FluxResult<LawType::Euler,nDim,Real> RoeUnprecWS::flux( const Species<LawType::Euler,Real>& species,
                                                           const geom::Surface<nDim,Real>&        face,
                                                           const StateT&                            sl,
                                                           const StateT&                            sr ) const
  {
      using FluxRes = FluxResult<LawType::Euler,
                                 nDim,
                                 Real>;

      using ConsDel = VariableDelta<LawType::Euler,
                                    nDim,
                                    BasisType<LawType::Euler>::Conserved,
                                    Real>;

      using ConsSet = VariableSet<LawType::Euler,
                                  nDim,
                                  BasisType<LawType::Euler>::Conserved,
                                  Real>;

      const StateT ravg = roeAverage( species, sl, sr );

   // central component
      const FluxRes central = CentralFlux<LawType::Euler>::flux( species, face, sl, sr );

   // diffusive component

      // interface jumps
      const Real dp = sr.pressure()
                     -sl.pressure();

      const Real du = projectedVelocity( face, sr )
                     -projectedVelocity( face, sl );

      // average quantities
      const Real una = projectedVelocity( face, ravg );
//    const Real ua  = sqrt( ravg.velocity2() );

      const Real ra = ravg.density();
      const Real aa = sqrt( ravg.speedOfSound2() );
//    const Real aa1= 1./aa;

//    const Real amu = aa - fabs( una );
//    const Real mna = una*aa1;
//    const Real ma  =  ua*aa1;

      const Real minf = species.minf;
      const Real uinf = minf*aa;

      // upwinding
      const ConsDel dqc =  state2Set<ConsSet>( species, sr )
                          -state2Set<ConsSet>( species, sl );

      const Real nu_p = (1./(ra*uinf*uinf))*pow( minf, 1 );
      const Real nu_u =                0.20*pow( minf, 0 );
      const Real mu_p =                0.00*pow( minf, 1 );
      const Real mu_u =           (ra*uinf)*pow( minf, 1 );

      // velocity perturbation deltaU*(rho,rhou,rhov,rhoh)
//    const Real deltaU =  amu*dp*aa1*aa1/ra
//                       +  ma*du;

      const Real deltaU =  nu_p*dp
                         + nu_u*du;

      const ConsDel ruvh = [&]() -> ConsDel
     {
         ConsDel qc;
         for( int i=0; i<nDim; ++i ){ qc[i]=ra*ravg.velocity(i); }
         qc[nDim  ]=ra;
         qc[nDim+1]=ra*ravg.specificTotalEnthalpy();
         return qc;
     }();

      // pressure perturbation deltaP*(0,nx,ny,U)
//    const Real deltaP =   ma*dp
//                       + amu*du*ra;

      const Real deltaP =  mu_p*dp
                         + mu_u*du;

      const ConsDel onu = [&]() -> ConsDel
     {
         ConsDel qc;
         for( int i=0; i<nDim; ++i ){ qc[i]=face.metric[0][i]; }
         qc[nDim  ]=0;
         qc[nDim+1]=una;
         return qc;
     }();

   // assemble

      const FluxRes diffusive{ 0.5*( fabs(una)*dqc
                                    +   deltaU*ruvh
                                    +   deltaP*onu )*face.area,
                               spectralRadius( face, ravg )*face.area};

      return FluxRes{ central.flux - diffusive.flux,
                      fmax( central.lambda, diffusive.lambda )};
  }
