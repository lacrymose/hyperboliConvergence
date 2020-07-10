
namespace geom
{
/*
 *  return identity metric
 */
   template<int nDim, floating_point Real>
   Metric<nDim,Real> metricI()
  {
      Metric<nDim,Real> identity{};
      for( int i=0; i<nDim; i++ ){ identity[i][i]=1.; }
      return identity;
  }

/*
 * return transpose metric
 *    metrics are unary matrices, so transpose is inverse
 */
   template<int nDim, floating_point Real>
   Metric<nDim,Real> transpose( const Metric<nDim,Real>& m )
  {
      Metric<nDim,Real> t;
      for( int i=0; i<nDim; i++ )
     {
         for( int j=0; j<nDim; j++ )
        {
            t[i][j] = m[j][i];
        }
     }
      return t;
  }
}
