
   template<floating_point Real>
   Species<LawType::Euler,Real> get_air_species()
  {
      Species<LawType::Euler,Real> air;

      air.gamma = 1.4;
      air.minf  = 0.0001;
      air.lref  = 1.;
      air.pr    = 0.7;
      air.nu    = 1.81e-5;
      air.dt    = std::numeric_limits<Real>::max();
      air.R     = 287.058;

      air.gamma1 = 1. / ( air.gamma - 1. );

      return air;
  }
