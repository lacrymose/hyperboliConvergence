
namespace Limiters
{
   template<typename Limiter>
   inline void ScalarFaceLimiter<Limiter>::operator()( Types::Real  dl,
                                                       Types::Real  dc,
                                                       Types::Real  dr,
                                                       Types::Real& dl_lim,
                                                       Types::Real& dr_lim ) const
  {
      dl_lim = limiter( dc, dl );
      dr_lim = limiter( dc, dr );
      return;
  }


   template<typename Limiter>
   template<typename DeltaType>
   inline void VectorFaceLimiter<Limiter>::operator()( const IdealGas2D::VariableDelta<DeltaType>& dl,
                                                       const IdealGas2D::VariableDelta<DeltaType>& dc,
                                                       const IdealGas2D::VariableDelta<DeltaType>& dr,
                                                             IdealGas2D::VariableDelta<DeltaType>& dl_lim,
                                                             IdealGas2D::VariableDelta<DeltaType>& dr_lim ) const
  {
      for( int i=0; i<4; i++ )
     {
         limiter( dl[i],dc[i],dr[i], dl_lim[i],dr_lim[i] );
     }
      return;
  }

   inline Types::Real NoLimit1::operator()( Types::Real a, Types::Real b ) const
  {
      return 0;
  }

   inline Types::Real NoLimit2::operator()( Types::Real a, Types::Real b ) const
  {
      return 0.5*(a+b);
  }

   inline Types::Real NoLimit3::operator()( Types::Real a, Types::Real b ) const
  {
      return 0.33333333333333333*(2.*a + b );
  }

   inline Types::Real VanAlbada2::operator()( Types::Real a, Types::Real b ) const
  {
      Types::Real z;

      z= ( b*(a*a) + a*(b*b) )/( a*a + b*b + Types::EPS );

      z*= a*b>0 ? 1:0;

      return z;
  };

   inline Types::Real MinMod2::operator()( Types::Real a, Types::Real b ) const
  {
      Types::Real       z;
      int           gt,lt;

      gt = a>0 ? 1:0;
      lt = a<0 ? 1:0;

      z =  gt*fmin( a, fmax( b, 0. ) )
         + lt*fmax( a, fmin( b, 0. ) );

      return z;
  };

   inline Types::Real VanLeer2::operator()( Types::Real a, Types::Real b ) const
  {
      return (b+fabs(b))/(a+fabs(b)+Types::EPS);
  }

   inline Types::Real Cada3::operator()( Types::Real a, Types::Real b ) const
  {

      Types::Real t =  0.0;
      Types::Real gam= 1.5;
      Types::Real c1= -2./7.;
      Types::Real c2=  0.4;
      Types::Real c3=  3*gam-2;
      Types::Real u,v,w,p,q;

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
