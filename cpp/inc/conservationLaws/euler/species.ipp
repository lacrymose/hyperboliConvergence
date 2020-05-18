
   inline Species<LawType::Euler> get_air_species()
  {
      Species<LawType::Euler> air;

      air.gamma = 1.4;
      air.minf  = 1.0;
      air.pr    = 0.7;
      air.nu    = 1.81e-5;
      air.dt    = std::numeric_limits<Types::Real>::max();
      air.R     = 2.87058;

      return air;
  }
