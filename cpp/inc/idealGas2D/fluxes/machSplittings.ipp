
namespace IdealGas2D
{
// mach number splittings for Van-Leer and AUSM type fluxes

   namespace MachSplitting
  {
      inline Types::Real M1( const int sign, const Types::Real M )
     {
         return 0.5*( M + sign*fabs(M) );
     }

      inline Types::Real M2( const int sign, const Types::Real M )
     {
         return sign*0.25*( M+sign )*( M+sign );
     }

      inline Types::Real M4( const int sign, const Types::Real M, const Types::Real beta )
     {
         return M2(sign,M)*( 1. - sign*16.*beta*M2(-sign,M) );
     }

      inline Types::Real P1( const int sign, const Types::Real M )
     {
         return (M1(sign,M)+Types::EPS)/(fabs(M)+Types::EPS);
     }

      inline Types::Real P3( const int sign, const Types::Real M )
     {
         return M2(sign,M)*( sign*2. - M );
     }

      inline Types::Real P5( const int sign, const Types::Real M, const Types::Real alpha )
     {
         return P3(sign,M) - sign*16.*alpha*M*M2(sign,M)*M2(-sign,M);
     }
  }
}
