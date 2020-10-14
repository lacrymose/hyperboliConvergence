
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

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
      std::array<nDim*3,Real> x;

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
   template<int nDim,
            floating_point Real>
   struct XMetricFactor
  {
      std::array<nDim*3,Real> x;

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
      static constexpr int ndim=dim_of_v<VarSet>;
      static constexpr int nvar=nVar<Law,ndim>;
      using VarSet = VarSet;
      using VarDel = vardelta_t<VarSet>;
      using Real = fptype_of_t<VarSet>;

      std::array<ndim*nvar,Real> q;

      int idx( const int i, const int j ){ return nvar*i+j; }

      const Real& operator( const int i, const int j ) const { return q[idx(i,j)]; }
            Real& operator( const int i, const int j )       { return q[idx(i,j)]; }

      QMetric& operator+=( const QMetric& other );
      QMetric& operator-=( const QMetric& other );
  };

   template<ImplementedVarSet VarSet>
   QMetric<VarSet> operator+( const QMetric<VarSet>& lhs,
                              const QMetric<VarSet>& rhs );

   template<ImplementedVarSet VarSet>
   QMetric<VarSet> operator-( const QMetric<VarSet>& lhs,
                              const QMetric<VarSet>& rhs );

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
   XMetric<nDim,Real> xmetric( const geom::Direction<nDim,Real>& d0 );
                     

/*
 *  calculation of Q*X^T metric for least squares calculation
 */
   template<int                 nDim,
            ImplementedVarSet VarSet,
            floating_point      Real>
      requires   nDim==dim_of_v<VarSet>
              && std::is_same_v<Real,fptype_of_t<VarSet>>
   QMetric<VarSet> qmetric( const geom::Point<nDim,Real>& x0,
                            const geom::Point<nDim,Real>& x1,
                            const VarSet&                 q0,
                            const VarSet&                 q1 )
  {
      return qmetric( x1-x0, q1-q0 );
  }

   template<int                   nDim,
            ImplementedVarDelta VarDel,
            floating_point        Real>
      requires   nDim==dim_of_v<VarDel>
              && std::is_same_v<Real,fptype_of_t<VarDel>>
   QMetric<VarSet> qmetric( const geom::Direction<nDim,Real>& d0,
                            const VarDel&                    dq0 );

}
