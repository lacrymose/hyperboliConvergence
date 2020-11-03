
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>
# include <utils/utils.h>

# include <utility>
# include <array>
# include <cmath>

namespace lsq
{

/*
 *  X*X^T metric for least squares calculation
 */
   template<int            nDim,
            floating_point Real>
   struct XMetric
  {
      std::array<utils::triangular_number(nDim),Real> x;

      const Real& operator( const int i ) const { return x[i]; }
            Real& operator( const int i )       { return x[i]; }

      XMetric& operator+=( const XMetric& other );
      XMetric& operator-=( const XMetric& other );
  };

   template<int            nDim,
            floating_point Real>
   XMetric<nDim,Real> operator+( const XMetric<nDim,Real>& lhs,
                                 const XMetric<nDim,Real>& rhs );

   template<int            nDim,
            floating_point Real>
   XMetric<nDim,Real> operator-( const XMetric<nDim,Real>& lhs,
                                 const XMetric<nDim,Real>& rhs );

/*
 *  factorised X*X^T metric for least squares calculation
 */
   template<int            nDim,
            floating_point Real>
   struct XMetricFactor
  {
      std::array<utils::triangular_number(nDim),Real> x;

      const Real& operator( const int i ) const { return x[i]; }
            Real& operator( const int i )       { return x[i]; }
  };

/*
 *  Q*X^T metric for least squares calculation
 */
   template<ImplementedVarSet VarSet>
   struct QMetric
  {
      static constexpr LawType Law = law_of_v<VarSet>;
      static constexpr int nDim=dim_of_v<VarSet>;
      static constexpr int nVar=nVar<Law,ndim>;
      using VarSet = varset_t<  VarSet>;
      using VarDel = vardelta_t<VarSet>;
      using Real = fptype_of_t< VarSet>;

      std::array<nDim*nVar,Real> q;

      int idx( const int i, const int j ){ return nVar*i+j; }

      const Real& operator( const int i, const int j ) const { return q[idx(i,j)]; }
            Real& operator( const int i, const int j )       { return q[idx(i,j)]; }

      QMetric& operator+=( const QMetric& other );
      QMetric& operator-=( const QMetric& other );
  };

   template<ImplementedVarSet VarSet>
   QMetric<VarT> operator+( const QMetric<VarSet>& lhs,
                            const QMetric<VarSet>& rhs );

   template<ImplementedVarSet VarSet>
   QMetric<VarT> operator-( const QMetric<VarSet>& lhs,
                            const QMetric<VarSet>& rhs );

/*
 * Convenience typedefs
 */
   // XMetric
   template<floating_point Real>
   using XMetric1 = XMetric<1,Real>;

   template<floating_point Real>
   using XMetric2 = XMetric<2,Real>;

   template<floating_point Real>
   using XMetric3 = XMetric<3,Real>;

   // XMetricFactor
   template<floating_point Real>
   using XMetricFactor1 = XMetricFactor<1,Real>;

   template<floating_point Real>
   using XMetricFactor2 = XMetricFactor<2,Real>;

   template<floating_point Real>
   using XMetricFactor3 = XMetricFactor<3,Real>;

/*
 *  calculation of X*X^T metric for least squares calculation
 */
   template<int            nDim,
            floating_point Real>
   XMetric<nDim,Real> xmetric( const geom::Point<nDim,Real>& x0,
                               const geom::Point<nDim,Real>& x1 )
  {
      return xmetric( x1-x0 );
  }

   template<int            nDim,
            floating_point Real>
   XMetric<nDim,Real> xmetric( const geom::Direction<nDim,Real>& dx );
                                                    
// template<floating_point Real>
// XMetric1<Real> xmetric( const geom::Direction1<Real>& dx );
//                                                
// template<floating_point Real>                  
// XMetric2<Real> xmetric( const geom::Direction2<Real>& dx );
//                                                
// template<floating_point Real>                  
// XMetric3<Real> xmetric( const geom::Direction3<Real>& dx );

/*
 * factors X*X^T metric
 */
   template<floating_point Real>
   XMetricFactor1<Real> factor( const XMetric1<Real>& xm );

   template<floating_point Real>
   XMetricFactor2<Real> factor( const XMetric2<Real>& xm );

   template<floating_point Real>
   XMetricFactor3<Real> factor( const XMetric3<Real>& xm );


/*
 *  calculation of Q*X^T metric for least squares calculation
 */
   template<int                 nDim,
            ImplementedVarSet VarSet,
            floating_point      Real>
      requires ConsistentTypes<law_of_v<VarSet>,
                               nDim,
                               Real,
                               VarSet>
   QMetric<VarSet> qmetric( const geom::Point<nDim,Real>& x0,
                            const geom::Point<nDim,Real>& x1,
                            const VarSet&                 q0,
                            const VarSet&                 q1 )
  {
      return qmetric( x1-x0,
                      q1-q0 );
  }

   template<int                   nDim,
            ImplementedVarDelta VarDel,
            floating_point        Real>
      requires ConsistentTypes<law_of_v<VarDel>,
                               nDim,
                               Real,
                               VarDel>
   QMetric<varset_t<VarDel>> qmetric( const geom::Direction<nDim,Real>& dx,
                                      const VarDel&                     dq );

// template<ImplementedVarDelta VarDel,
//          floating_point        Real>
//    requires ConsistentTypes<law_of_v<VarDel>,
//                             1,
//                             Real,
//                             VarDel>
// QMetricVarDel> qmetric( const geom::Direction1<Real>& dx,
//                         const VarDel&                 dq );

// template<ImplementedVarDelta VarDel,
//          floating_point        Real>
//    requires ConsistentTypes<law_of_v<VarDel>,
//                             2,
//                             Real,
//                             VarDel>
// QMetric<VarDel> qmetric( const geom::Direction2<Real>& dx,
//                          const VarDel&                 dq );

// template<ImplementedVarDelta VarDel,
//          floating_point        Real>
//    requires ConsistentTypes<law_of_v<VarDel>,
//                             3,
//                             Real,
//                             VarDel>
// QMetric<VarDel> qmetric( const geom::Direction3<Real>& dx,
//                          const VarDel&                 dq );

/*
 * solve for the gradient in a given direction
 */
   template<ImplementedVarSet VarSet,
            int                 nDim,
            floating_point      Real>
      requires ConsistentTypes<law_of_v<VarSet>,
                               nDim,
                               Real,
                               VarSet>
   vardelta_t<VarSet> solve( const XMetric<nDim,Real>&          xm,
                             const QMetric<VarSet>&             qm,
                             const geom::Direction<nDim,Real>& dir )
  {
      return solve( factor( xm ), qm, dir );
  }

   template<ImplementedVarSet VarSet,
            floating_point      Real>
      requires ConsistentTypes<law_of_v<VarSet>,
                               1,
                               Real,
                               VarSet>
   vardelta_t<VarSet> solve( const XMetricFactor1<Real>&    xf,
                             const QMetric<VarSet>&         qm,
                             const geom::Direction1<Real>& dir );

/*
 * solve for the gradient in a given direction
 */
   template<ImplementedVarSet VarSet,
            floating_point      Real>
      requires ConsistentTypes<law_of_v<VarSet>,
                               2,
                               Real,
                               VarSet>
   vardelta_t<VarSet> solve( const XMetricFactor2<Real>&    xf,
                             const QMetric<VarSet>&         qm,
                             const geom::Direction2<Real>& dir );

}

# include <lsq/xmetric.ipp>
# include <lsq/qmetric.ipp>
# include <lsq/factor_solve.ipp>
