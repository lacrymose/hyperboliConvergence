
# include <limits>

// mach number splittings for Van-Leer and AUSM type fluxes

   template<floating_point Real>
   Real MachSplit_M1( const int sign, const Real M )
  {
      return 0.5*( M + sign*fabs(M) );
  }

   template<floating_point Real>
   Real MachSplit_M2( const int sign, const Real M )
  {
      return sign*0.25*( M+sign )*( M+sign );
  }

   template<floating_point Real>
   Real MachSplit_M4( const int sign, const Real M, const Real beta )
  {
      return MachSplit_M2(sign,M)*( 1. - sign*16.*beta*MachSplit_M2(-sign,M) );
  }

   template<floating_point Real>
   Real MachSplit_P1( const int sign, const Real M )
  {
      constexpr Real eps=std::numeric_limits<Real>::max();
      return (MachSplit_M1(sign,M)+eps)/(fabs(M)+eps);
  }

   template<floating_point Real>
   Real MachSplit_P3( const int sign, const Real M )
  {
      return MachSplit_M2(sign,M)*( sign*2. - M );
  }

   template<floating_point Real>
   Real MachSplit_P5( const int sign, const Real M, const Real alpha )
  {
      return MachSplit_P3(sign,M) - sign*16.*alpha*M*MachSplit_M2(sign,M)*MachSplit_M2(-sign,M);
  }
