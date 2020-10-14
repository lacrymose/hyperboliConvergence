
# pragma once

# include <utils/concepts.h>
# include <utils/maths/misc.h>

# include <array>

namespace fem
{
// 1D functions

   // shape functions
   template<floating_point Real>
   Real n10( Real x ){ return 0.5*( 1.-x ); }

   template<floating_point Real>
   Real n11( Real x ){ return 0.5*( 1.+x ); }

   template<floating_point Real>
   utils::matrix_1<Real,2> n1( Real x ){ return {n10(x),
                                                 n11(x)}; }

   // derivatives
   template<floating_point Real>
   Real dn10( Real x ){ return -0.5; }

   template<floating_point Real>
   Real dn11( Real x ){ return  0.5; }

   template<floating_point Real>
   utils::matrix_1<Real,2> dn1( Real x ){ return {dn10(x),
                                                  dn11(x)}; }

// 2D functions

   // shape functions
   template< floating_point Real>
   Real n200( utils::matrix_1<Real,2> x ){ return n10(x[0])*n10(x[1]); }

   template< floating_point Real>
   Real n210( utils::matrix_1<Real,2> x ){ return n11(x[0])*n10(x[1]); }

   template< floating_point Real>
   Real n201( utils::matrix_1<Real,2> x ){ return n10(x[0])*n11(x[1]); }

   template< floating_point Real>
   Real n211( utils::matrix_1<Real,2> x ){ return n11(x[0])*n11(x[1]); }

   template< floating_point Real>
   utils::matrix_1<Real,4> n2( utils::matrix_1<Real,2> x ){ return{n200(x),
                                                                   n210(x),
                                                                   n201(x),
                                                                   n211(x)}; }

   // derivatives
   template< floating_point Real>
   utils::matrix_1<Real,2> dn200( utils::matrix_1<Real,2> x ){ return {dn10(x[0])* n10(x[1]),
                                                                        n10(x[0])*dn10(x[1])}; }

   template< floating_point Real>
   utils::matrix_1<Real,2> dn210( utils::matrix_1<Real,2> x ){ return {dn11(x[0])* n10(x[1]),
                                                                        n11(x[0])*dn10(x[1])}; }

   template< floating_point Real>
   utils::matrix_1<Real,2> dn201( utils::matrix_1<Real,2> x ){ return {dn10(x[0])* n11(x[1]),
                                                                        n10(x[0])*dn11(x[1])}; }

   template< floating_point Real>
   utils::matrix_1<Real,2> dn211( utils::matrix_1<Real,2> x ){ return {dn11(x[0])* n11(x[1]),
                                                                        n11(x[0])*dn11(x[1])}; }

   template< floating_point Real>
   utils::matrix_2<Real,4,2> dn2( utils::matrix_1<Real,2> x ){ return {dn200(x),
                                                                       dn210(x),
                                                                       dn201(x),
                                                                       dn211(x)}; }
}
