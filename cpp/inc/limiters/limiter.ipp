
# include <limits>

# include <cmath>

namespace Limiters
{
   template<floating_point Real>
   Real NoLimit1::operator()( const Real a, const Real b ) const
  {
      return 0;
  }

   template<floating_point Real>
   Real NoLimit2::operator()( const Real a, const Real b ) const
  {
      return 0.5*(a+b);
  }

   template<floating_point Real>
   Real NoLimit3::operator()( const Real a, const Real b ) const
  {
      return 0.33333333333333333*(2.*a + b );
  }

   template<floating_point Real>
   Real VanAlbada2::operator()( const Real a, const Real b ) const
  {
      Real z;

      const Real eps = std::numeric_limits<Real>::min();

   // do not evaluate division if unnecessary
      z = a*b<0 ?  0
                 : ( b*(a*a) + a*(b*b) )
                  /(    a*a  +    b*b + eps );

      return z;
  }

   template<floating_point Real>
   Real MinMod2::operator()( const Real a, const Real b ) const
  {
      Real       z;

   // do not evaluate extra min/max if unnecessary
      z = a>0 ?  fmin( a, fmax( b, 0. ) )
               : fmax( a, fmin( b, 0. ) );

      return z;
  }

   template<floating_point Real>
   Real VanLeer2::operator()( const Real a, const Real b ) const
  {
      const Real eps = std::numeric_limits<Real>::min();
      return ( b + fabs(b) )
            /( a + fabs(b) + eps );
  }

   template<floating_point Real>
   Real Cada3::operator()( const Real a, const Real b ) const
  {

      Real t =  0.0;
      Real gam= 1.5;
      Real c1= -2./7.;
      Real c2=  0.4;
      Real c3=  3*gam-2;
      Real u,v,w,p,q;

      if( a > 0 )
     {
         u= a;
         v= b;
     }
      else
     {
         u=-a;
         v=-b;
     }
      p= (v+ 2.*u)/3.;
      q= fabs(a)+ fabs(b); 
      if( q > t )
     {

         if( v < -2*u )
        {
            w=0;
        }
         else
        {
            if( v < c1*u )
           {
              w= (v+2.*u)/3.;
           }
            else
           {
               if( v < 0 )
              {
                  w= -2.*v;
              }
               else
              {
                  if( v < c2*u )
                 {
                      w= 2.*v;
                 }
                  else
                 {
                     if( v < c3*u )
                    {
                        w= (v+2.*u)/3.;
                    }
                     else
                    {
                         w=gam*u;
                    }
                 }
              }
           }
        }
     }
      else
     {
         w= p;
     }
      if( a <= 0 ){ w= -w; };

      return w;
  }
}
