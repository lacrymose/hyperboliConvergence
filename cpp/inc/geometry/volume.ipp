
# include <fem/fem.h>
# include <utils/maths/misc.h>

# include <cassert>

namespace geom
{
   template<floating_point Real>
   Volume<1,Real> volume( const Point<1,Real>& p0, const Point<1,Real>& p1 )
  {
      const Direction<1,Real> dir=p1-p0;
      return {.volume=length(dir),
              .centre=p0+0.5*dir};
  }

   // corners of standard quad
   //
   //   2 ---- 3
   //   |      |
   //   |      |
   //   0 ---- 1
   //

   template<floating_point Real>
   Volume<2,Real> volume( const Point<2,Real>& p0, const Point<2,Real>& p1,
                          const Point<2,Real>& p2, const Point<2,Real>& p3 )
  {
   // shape function and derivative
      constexpr Real zero=0.;
      const utils::matrix_1<Real,4>    u = fem::n2<Real>(  {zero,zero} );
      const utils::matrix_2<Real,4,2> du = fem::dn2<Real>( {zero,zero} );

   // centre of mass
      Point<2,Real> com{};
      for( size_t i=0; i<2; ++i )
     {
         com[i]+= u[0]*p0[i];
         com[i]+= u[1]*p1[i];
         com[i]+= u[2]*p2[i];
         com[i]+= u[3]*p3[i];
     }

   // area of quadrilateral
      std::array<Direction<2,Real>,2> tangent{};
      for( size_t i=0; i<2; ++i )
     {
         for( size_t j=0; j<2; ++j )
        {
            tangent[i][j]+=du[0][i]*p0[j];
            tangent[i][j]+=du[1][i]*p1[j];
            tangent[i][j]+=du[2][i]*p2[j];
            tangent[i][j]+=du[3][i]*p3[j];
        }
     }

      const Real vol = std::fabs( cross( tangent[0],
                                         tangent[1] ) );

      assert( vol>0. );
      assert( !std::isnan(vol) );
      assert(  std::isfinite(vol) );

      return {.volume=vol,
              .centre=com};
  }
}
