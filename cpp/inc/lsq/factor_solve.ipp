
namespace lsq
{
/*
 * factors X*X^T metric
 */
   template<floating_point Real>
   XMetricFactor1<Real> factor( const XMetric1<Real>& xm )
  {
      return XMetricFactor1<Real>{{1./xm(0)}};
  }

   template<floating_point Real>
   XMetricFactor2<Real> factor( const XMetric2<Real>& xm )
  {
      XMetricFactor2<Real> xf;
      const Real d = 1./( xm(0)*xm(2) - xm(1)*xm(1) );
      xf(0) =  d*xm(2);
      xf(1) = -d*xm(1);
      xf(2) =  d*xm(0);
      return xf;
  }

/*
 * solve for the gradient in a given direction
 */
   template<ImplementedVarSet VarSet,
            floating_point      Real>
      requires ConsistentTypes<law_of_v<VarSet>,
                               1,
                               Real,
                               VarSet>
   vardelta_t<VarSet> solve( const XMetricFactor1<Real>&    xf,
                             const QMetric<VarSet>&         qm,
                             const geom::Direction1<Real>& dir )
  {
      using VarDel = vardelta_t<VarSet>;
      const auto d = geom::norm( dir );
      VarDel dq;
      for( unsigned int i=0; i<VarDel::N; ++i )
     {
         dq[i] = qm(i,0)*xf(0)*dir[0];
     }
      return dq;
  }

   template<ImplementedVarSet VarSet,
            floating_point      Real>
      requires ConsistentTypes<law_of_v<VarSet>,
                               2,
                               Real,
                               VarSet>
   vardelta_t<VarSet> solve( const XMetricFactor2<Real>&    xf,
                             const QMetric<VarSet>&         qm,
                             const geom::Direction2<Real>& dir )
  {
      using VarDel = vardelta_t<VarSet>;
      const auto d = geom::norm( dir

      VarDel dq;
   // xf*qf is gradient (vector quantity)
   // xf*qf*d is gradient in direction d (scalar)
      for( unsigned int i=0; i<VarDel::N; ++i )
     {
         dq[i] =  ( xf(0)*qm(i,0) + xf(1)*qm(i,1) )*dir[0]
                + ( xf(1)*qm(i,0) + xf(2)*qm(i,1) )*dir[1];
     }
      return dq;
  }

   template<ImplementedVarSet VarSet,
            int                 nDim,
            floating_point      Real>
      requires ConsistentTypes<law_of_v<VarSet>,
                               nDim,
                               Real,
                               VarSet>
   std::pair<VarDel,VarDel> bias_solves( const geom::Direction<nDim,Real>& dx_c,
                                         const VarDel&                     dq_c,
                                         const XMetric<nDim,Real>&        dxm_l,
                                         const XMetric<nDim,Real>&        dxm_r,
                                         const QMetric<VarSet>&           dqm_l,
                                         const QMetric<VarSet>&           dqm_r,
                                         const geom::Direction<nDim,Real>&  dir )
  {
      // biased left/right spatial metrics
         const auto dxm_c  = lsq::xmetric(dx_c);
         const auto dxmb_l = dxm_l - dxm_c;
         const auto dxmb_r = dxm_r - dxm_c;

      // biased left/right solution metrics
         const auto dqm_c  = lsq::qmetric(dx_c, dq_c);
         const auto dqmb_l = dqm_l - dqm_c;
         const auto dqmb_r = dqm_r - dqm_c;

      // l/r biased gradients
      //    long 1-liner to make use of RVO
         return {solve( dxmb_l, dqmb_l, dir ),
                 solve( dxmb_r, dqmb_r, dir )};
  }
}
