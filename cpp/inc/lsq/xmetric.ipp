
namespace lsq
{

/*
 * In-place arithmetic is element-wise
 */
   template<int            nDim,
            floating_point Real>
   XMetric<nDim,Real>& XMetric<nDim,Real>::operator+=( const XMetric<nDim,Real>& other )
  {
      for( unsigned int i=0; i<x.size(); ++i )
     {
         x[i]+=other.x[i];
     }
      return *this;
  }

   template<int            nDim,
            floating_point Real>
   XMetric<nDim,Real>& XMetric<nDim,Real>::operator-=( const XMetric<nDim,Real>& other )
  {
      for( unsigned int i=0; i<x.size(); ++i )
     {
         x[i]-=other.x[i];
     }
      return *this;
  }

   template<int            nDim,
            floating_point Real>
   XMetric<nDim,Real> operator+( const XMetric<nDim,Real>& lhs,
                                 const XMetric<nDim,Real>& rhs )
  {
      XMetric<nDim,Real> result(lhs);
      result+=rhs;
      return result;
  }

   template<int            nDim,
            floating_point Real>
   XMetric<nDim,Real> operator-( const XMetric<nDim,Real>& lhs,
                                 const XMetric<nDim,Real>& rhs )
  {
      XMetric<nDim,Real> result(lhs);
      result-=rhs;
      return result;
  }

/*
 *  calculation of X*X^T metric for least squares calculation
 */
   template<int            nDim,
            floating_point Real>
   XMetric<nDim,Real> xmetric( const geom::Direction<nDim,Real>& dx )
  {
      XMetric<nDim,Real> xm;

      unsigned int k=0;
      for( unsigned int i=0; i<nDim; ++i )
     {
         for( unsigned int j=i; j<nDim; ++j )
        {
            xm(k) = dx[i]*dx[j];
            ++k;
        }
     }
      return xm;
  }
}
