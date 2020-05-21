
   template<floating_point Real>
   Species<LawType::Euler,Real> get_air_species()
  {
      Species<LawType::Euler,Real> air;

      air.gamma = 1.4;
      air.minf  = 1.0;
      air.pr    = 0.7;
      air.nu    = 1.81e-5;
      air.dt    = std::numeric_limits<Real>::max();
      air.R     = 2.87058;

      air.gamma1 = 1. / ( air.gamma - 1. );

      return air;
  }
