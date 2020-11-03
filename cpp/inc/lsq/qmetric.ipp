
namespace lsq
{
/*
 * In-place arithmetic is element-wise
 */
   template<ImplementedVarSet VarSet>
   QMetric<VarSet>& QMetric<VarSet>::operator+=( const QMetric<VarSet>& other )
  {
      for( unsigned int i=0; i<q.size(); ++i )
     {
         q[i]+=other.q[i];
     }
      return *this;
  }

   template<ImplementedVarSet VarSet>
   QMetric<VarSet>& QMetric<VarSet>::operator-=( const QMetric<VarSet>& other )
  {
      for( unsigned int i=0; i<q.size(); ++i )
     {
         q[i]-=other.q[i];
     }
      return *this;
  }

   template<ImplementedVarSet VarSet>
   QMetric<VarSet> operator+( const QMetric<VarSet>& lhs,
                              const QMetric<VarSet>& rhs )
  {
      QMetric<VarSet> result(lhs);
      result+=rhs;
      return result;
  }

   template<ImplementedVarSet VarSet>
   QMetric<VarSet> operator-( const QMetric<VarSet>& lhs,
                              const QMetric<VarSet>& rhs )
  {
      QMetric<VarSet> result(lhs);
      result-=rhs;
      return result;
  }

/*
 *  calculation of Q*X^T metric for least squares calculation
 */
   template<int                   nDim,
            ImplementedVarDelta VarDel,
            floating_point        Real>
      requires ConsistentTypes<law_of_v<VarDel>,
                               nDim,
                               Real,
                               VarDel>
   QMetric<varset_t<VarDel>> qmetric( const geom::Direction<nDim,Real>& dx,
                                      const VarDel&                     dq )
  {
      using QMet = QMetric<varset_t<VarDel>>;

      QMet qm;

      for( unsigned int i=0; i<QMet::nVar; ++i )
     {
         for( unsigned int j=0; j<QMet::nDim; ++j )
        {
            qm(i,j) = dq[i]*dx[j];
        }
     }
      return qm;
  }
}
