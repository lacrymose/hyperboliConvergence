

   template<LawType Law>
   template<int nDim, floating_point Real>
   FluxResult<Law,nDim,Real> RoeFlux<Law>::flux( const Species<Law,Real>&     species,
                                                 const geom::Surface<nDim,Real>& face,
                                                 const State<Law,nDim,Real>&       sl,
                                                 const State<Law,nDim,Real>&       sr )
  {

   // types needed for upwind diffusion calculation
      using ConsSet   = VariableSet<  Law,nDim,BasisType<Law>::Conserved,     Real>;
      using ConsDelta = VariableDelta<Law,nDim,BasisType<Law>::Conserved,     Real>;
      using CharDelta = VariableDelta<Law,nDim,BasisType<Law>::Characteristic,Real>;
      using FluxRes   = FluxResult<Law,nDim,Real>;
      using StateT    = State<Law,nDim,Real>;

   // central flux
      FluxRes result = CentralFlux<Law>::flux( species, face, sl,sr );

   // matrix upwind diffusion

      // roe average state and eigenvalues
      const StateT ravg = roeAverage( species, sl, sr );

      const WaveSpeeds<Law,nDim,Real> lambda = entropyfix( wavespeeds( species, face, ravg ),
                                                           wavespeeds( species, face,  sl  ),
                                                           wavespeeds( species, face,  sr  ) );

   // face aligned conserved variables
      const ConsDelta dqc = rotateToFace( face, ( state2Set<ConsSet>( species, sr )
                                                 -state2Set<ConsSet>( species, sl ) ) );

   // scaled characteristic variables
      const CharDelta dqw = lambda*delta2Delta<CharDelta>( species, ravg, dqc );

   // roe dissipation in conserved variables
      const ConsDelta fd = rotateFromFace( face, delta2Delta<ConsDelta>( species, ravg, dqw ) );

      // spectral radius
      result.lambda = fmax( result.lambda,
                            spectralRadius( face, ravg )*face.area );

      result.flux-=0.5*fd;

      return result;
  }


