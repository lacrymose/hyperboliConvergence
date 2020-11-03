
# pragma once

# include <lsq/lsq.h>

# include <conservationLaws/base/base.h>

# include <mesh/mesh.h>
# include <geometry/geometry.h>

# include <parallalg/neighbour_algorithm.h>
# include <parallalg/array.h>
# include <parallalg/parallalg.h>

# include <functional>

/*
 * Calculate array of spatial metrics for least squares gradient calculation
 */
   template<par::execution_policy Policy,
            int                     nDim,
            floating_point          Real>
   auto xmetrics( const Policy                   policy,
                  const MeshCellArray<nDim,Real>& cells )
  {
      par::DualArray<lsq::XMetric<nDim,Real>,nDim> dxdx(cells.shape());
      xmetrics( policy, cells, dxdx );
      return dxdx;
  }

   template<par::execution_policy Policy,
            int                     nDim,
            floating_point          Real>
   void xmetrics( const Policy                                      policy,
                  const MeshCellArray<nDim,Real>&                    cells,
                        par::DualArray<lsq::XMetric<nDim,Real>,nDim>& dxdx )
  {
      using XMetric = lsq::XMetric<nDim,Real>;
      using Cell = typename Mesh<nDim,Real>::Cell;

      assert( cells.shape() == dxdx.shape() );

      par::fill( policy, dxdx, XMetric{0.} )

      const auto xmetric_calc = []( const Cell& c0,
                                    const Cell& c1 ) -> XMetric
     {
         return lsq::xmetric( c0.centre,
                              c1.centre );
     }

      par::neighbour_accumulation( policy,
                                   xmetric_calc,
                                   std::plus<XMetric>,
                                   std::plus<XMetric>,
                                   dxdx,
                                   cells );

      boundary_xmetrics( policy, cells, dxdx );

      return;
  }

/*
 * use boundary contributions to xmetrics to prevent rank deficient biased metrics
 */
   template<par::execution_policy Policy,
            floating_point          Real>
   void boundary_xmetrics( const Policy                                 policy,
                           const MeshCellArray<1,Real>&                  cells,
                                 par::DualArray<lsq::XMetric<1,Real>,1>&  dxdx )
  {
      using Cell = typename  MeshCellArray<1,Real>::ElemType;
      using XMetric = lsq::XMetric<1,Real>;

      const auto boundary_dxdx_calc = []( const Cell& c0,
                                                Cell& c1 ) -> XMetric
     {
         const auto displacement = cell0.centre - cell1.centre;
         const auto boundary_point = cell0.centre + displacement;

        return lsq::xmetric( cell0.centre,
                             boundary_point );
     };

   // left face
     {
         const par::DualIdx1 i0{0};
         const par::DualIdx1 i1{1};

         dxdx(i0) += boundary_dxdx_calc( cells(i0), cells(i1) );
     }
   // right face
     {
         const auto nc=cells.shape(0);
         const par::DualIdx1 i0{nc-1};
         const par::DualIdx1 i1{nc-2};

         dxdx(i0) += boundary_dxdx_calc( cells(i0), cells(i1) );
     }
      return;
  }

   template<par::execution_policy Policy,
            floating_point          Real>
   void boundary_xmetrics( const Policy                                 policy,
                           const MeshCellArray<2,Real>&                  cells,
                                 par::DualArray<lsq::XMetric<2,Real>,2>&  dxdx )
  {
      const auto ni = cells.shape(0);
      const auto nj = cells.shape(1);

      using Cell = typename  MeshCellArray<2,Real>::ElemType;
      using XMetric = lsq::XMetric<2,Real>;

      const auto boundary_dxdx_calc = []( const Cell& c0,
                                                Cell& c1 ) -> XMetric
     {
         const auto displacement = cell0.centre - cell1.centre;
         const auto boundary_point = cell0.centre + displacement;

        return lsq::xmetric( cell0.centre,
                             boundary_point );
     };

   // left/right faces
      for( size_t j=0; j<nj; ++j )
     {
      // left face
        {
            const par::DualIdx2 i0{0,j};
            const par::DualIdx2 i1{1,j};

            dxdx(i0) += boundary_dxdx_calc( cells(i0), cells(i1) );
        }
      // right face
        {
            const par::DualIdx2 i0{ni-1,j};
            const par::DualIdx2 i1{ni-2,j};

            dxdx(i0) += boundary_dxdx_calc( cells(i0), cells(i1) );
        }
     }
   // bottom/top faces
      for( size_t i=0; i<ni; ++i )
     {
      // bottom face
        {
            const par::DualIdx2 i0{i,0};
            const par::DualIdx2 i1{i,1};

            dxdx(i0) += boundary_dxdx_calc( cells(i0), cells(i1) );
        }
      // top face
        {
            const par::DualIdx2 i0{i,nj-1};
            const par::DualIdx2 i1{i,nj-2};

            dxdx(i0) += boundary_dxdx_calc( cells(i0), cells(i1) );
        }
     }
      return;
  }

/*
 * Calculate array of solution metrics for least squares gradient calculation
 */
   template<par::execution_policy  Policy,
            ImplementedVarSet   SolVarSet,
            int                      nDim,
            floating_point           Real>
   auto qmetrics( const Policy                        policy,
                  const MeshCellArray<nDim,Real>&      cells,
                        par::DualArray<SolVarSet,nDim>&   q0 )
  {
      par::DualArray<lsq::QMetric<SolVarSet>,nDim> dqdx(cells.shape());
      qmetrics( policy, cells, q0, dqdx );
      return dqdx;
  }

   template<par::execution_policy   Policy,
            ImplementedVarSet    SolVarSet,
            int                       nDim,
            floating_point            Real>
   void qmetrics( const Policy                                      policy,
                  const MeshCellArray<nDim,Real>&                    cells,
                  const par::DualArray<SolVarSet,nDim>&                 q0,
                        par::DualArray<lsq::QMetric<SolVarSet>,nDim>& dqdx )
  {
      using QMetric = lsq::QMetric<SolVarSet>;
      using Cell = typename MeshCellArray<nDim,Real>::ElemType;

      par::fill( policy, dqdx, QMetric{0.} )

      const auto qmetric_calc = []( const Cell&      c0,
                                    const Cell&      c1,
                                    const SolVarSet& q0,
                                    const SolVarSet& q1 ) ->QMetric
     {
         return lsq::qmetric( c0.centre,
                              c1.centre,
                              q0,
                              q1 );
     }

      par::neighbour_accumulation( policy,
                                   qmetric_calc,
                                   std::plus<QMetric>,
                                   std::plus<QMetric>,
                                   dqdx,
                                   cells, q0 );
      return;
  }

