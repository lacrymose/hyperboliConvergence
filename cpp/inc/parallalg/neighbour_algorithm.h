
# pragma once

# include <parallalg/array.h>
# include <parallalg/parallalg.h>

# include <omp.h>

namespace par
{
/*
 * run in serial if no policy provided
 */
   template<typename       EdgeFuncObj,
            typename    AccLeftFuncObj,
            typename   AccRightFuncObj,
            typename            ElemTd,
            typename            ElemT0,
            typename...         ElemTs,
            int                   NDIM,
            GridType               GTd,
            GridType               GT0,
            GridType...            GTs,
            ArraySizing            ASd,
            ArraySizing            AS0,
            ArraySizing...         ASs>
   void neighbour_accumulation(       EdgeFuncObj                 edge_func,
                                      AccLeftFuncObj              accl_func,
                                      AccRightFuncObj             accr_func,
                                      Array<ElemTd,NDIM,GTd,ASd>&       dst,
                                const Array<ElemT0,NDIM,GT0,AS0>&      src0,
                                const Array<ElemTs,NDIM,GTs,ASs>&...   srcs )
  {
      neighbour_accumulation( execution::seq,
                              edge_func, accl_func, accr_func,
                              dst, src0, srcs... );
      return;
  }

/*
 * apply a function to a 2-element stencil from a list of arrays
 */
   template<typename               FuncObj,
            int                       NDIM,
            typename                 ElemT,
            GridType                    GT,
            ArraySizing                 AS,
            typename... ArraysAndArguments>
   auto apply_stencil2(       FuncObj          stencil2_func,
                        const Idx<NDIM,GT>&             idxl,
                        const Idx<NDIM,GT>&             idxr,
                        const Array<ElemT,NDIM,GT,AS>&  src0,
                        const ArraysAndArguments...     args ) _PAR_ALWAYS_INLINE_
  {
      return apply_stencil2( stencil2_func,
                             idxl, idxr,
                             args...,
                             src0(idxl), src0(idxr) );
  }

   template<typename      FuncObj,
            int              NDIM,
            GridType           GT,
            typename... Arguments>
   auto apply_stencil2(       FuncObj   stencil2_func,
                        const Idx<NDIM,GT>&      idxl,
                        const Idx<NDIM,GT>&      idxr,
                        const Arguments...       args ) _PAR_ALWAYS_INLINE_
  {
      return stencil2_func( args... );
  }

/*
 * 1D all arrays are same grid type - serial execution
 */
   template<typename       EdgeFuncObj,
            typename    AccLeftFuncObj,
            typename   AccRightFuncObj,
            typename            ElemTd,
            typename            ElemT0,
            typename...         ElemTs,
            GridType                GT,
            ArraySizing            ASd,
            ArraySizing            AS0,
            ArraySizing...         ASs>
   void neighbour_accumulation(       execution::serial_policy,
                                      EdgeFuncObj             edge_func,
                                      AccLeftFuncObj          accl_func,
                                      AccRightFuncObj         accr_func,
                                      Array1<ElemTd,GT,ASd>&        dst,
                                const Array1<ElemT0,GT,AS0>&       src0,
                                const Array1<ElemTs,GT,ASs>&...    srcs )
  {
         assert(   (dst.shape() == src0.shape())
                && "par::neighbour_accumulation - arrays must be the same shape" );

      (( assert(   (dst.shape() == srcs.shape())
                && "par::neighbour_accumulation - arrays must be the same shape" ) ),... );

      const size_t ni = dst.shape(0);

      for( size_t i=0; i<ni-1; ++i )
     {
         const Idx1<GT> il{i  };
         const Idx1<GT> ir{i+1};

         const auto acc_v = apply_stencil2( edge_func,
                                            il,ir,
                                            src0, srcs... );

         dst(il) = accl_func( std::move(dst(il)),
                              acc_v );

         dst(ir) = accr_func( std::move(dst(ir)),
                              acc_v );
     }
  }

/*
 * 1D all arrays are same grid type - OpenMP execution
 */
   template<typename       EdgeFuncObj,
            typename    AccLeftFuncObj,
            typename   AccRightFuncObj,
            typename            ElemTd,
            typename            ElemT0,
            typename...         ElemTs,
            GridType                GT,
            ArraySizing            ASd,
            ArraySizing            AS0,
            ArraySizing...         ASs>
   void neighbour_accumulation(       execution::openmp_policy,
                                      EdgeFuncObj             edge_func,
                                      AccLeftFuncObj          accl_func,
                                      AccRightFuncObj         accr_func,
                                      Array1<ElemTd,GT,ASd>&        dst,
                                const Array1<ElemT0,GT,AS0>&       src0,
                                const Array1<ElemTs,GT,ASs>&...    srcs )
  {
         assert(   (dst.shape() == src0.shape())
                && "par::neighbour_accumulation - arrays must be the same shape" );

      (( assert(   (dst.shape() == srcs.shape())
                && "par::neighbour_accumulation - arrays must be the same shape" ) ),... );

      const size_t ni = dst.shape(0);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni-1; ++i )
     {
         const Idx1<GT> il{i  };
         const Idx1<GT> ir{i+1};

         const auto acc_v = apply_stencil2( edge_func,
                                            il,ir,
                                            src0, srcs... );

         dst(il) = accl_func( std::move(dst(il)),
                              acc_v );

         dst(ir) = accr_func( std::move(dst(ir)),
                              acc_v );
     }
  }

/*
 * 2D all arrays are same grid type - serial execution
 */
   template<typename       EdgeFuncObj,
            typename    AccLeftFuncObj,
            typename   AccRightFuncObj,
            typename            ElemTd,
            typename            ElemT0,
            typename...         ElemTs,
            GridType                GT,
            ArraySizing            ASd,
            ArraySizing            AS0,
            ArraySizing...         ASs>
   void neighbour_accumulation(       execution::serial_policy,
                                      EdgeFuncObj             edge_func,
                                      AccLeftFuncObj          accl_func,
                                      AccRightFuncObj         accr_func,
                                      Array2<ElemTd,GT,ASd>&        dst,
                                const Array2<ElemT0,GT,AS0>&       src0,
                                const Array2<ElemTs,GT,ASs>&...    srcs )
  {
         assert(   (dst.shape() == src0.shape())
                && "par::neighbour_accumulation - arrays must be the same shape" );

      (( assert(   (dst.shape() == srcs.shape())
                && "par::neighbour_accumulation - arrays must be the same shape" ) ),... );

      const size_t ni = dst.shape(0);
      const size_t nj = dst.shape(1);

      for( size_t j=0; j<nj; ++j )
     {
         for( size_t i=0; i<ni-1; ++i )
        {
            const Idx2<GT> ijl{i  ,j};
            const Idx2<GT> ijr{i+1,j};

            const auto acc_v = apply_stencil2( edge_func,
                                               ijl,ijr,
                                               src0, srcs... );

            dst(ijl) = accl_func( std::move(dst(ijl)),
                                  acc_v );

            dst(ijr) = accr_func( std::move(dst(ijr)),
                                  acc_v );
        }
     }

   // derivatives over i-normal faces
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj-1; ++j )
        {
            const Idx2<GT> ijl{i,j  };
            const Idx2<GT> ijr{i,j+1};

            const auto acc_v = apply_stencil2( edge_func,
                                               ijl,ijr,
                                               src0, srcs... );

            dst(ijl) = accl_func( std::move(dst(ijl)),
                                  acc_v );

            dst(ijr) = accr_func( std::move(dst(ijr)),
                                  acc_v );
        }
     }

      return;
  }

/*
 * 2D all arrays are same grid type - OpenMP execution
 */
   template<typename       EdgeFuncObj,
            typename    AccLeftFuncObj,
            typename   AccRightFuncObj,
            typename            ElemTd,
            typename            ElemT0,
            typename...         ElemTs,
            GridType                GT,
            ArraySizing            ASd,
            ArraySizing            AS0,
            ArraySizing...         ASs>
   void neighbour_accumulation(       execution::openmp_policy,
                                      EdgeFuncObj             edge_func,
                                      AccLeftFuncObj          accl_func,
                                      AccRightFuncObj         accr_func,
                                      Array2<ElemTd,GT,ASd>&        dst,
                                const Array2<ElemT0,GT,AS0>&       src0,
                                const Array2<ElemTs,GT,ASs>&...    srcs )
  {
         assert(   (dst.shape() == src0.shape())
                && "par::neighbour_accumulation - arrays must be the same shape" );

      (( assert(   (dst.shape() == srcs.shape())
                && "par::neighbour_accumulation - arrays must be the same shape" ) ),... );

      const size_t ni = dst.shape(0);
      const size_t nj = dst.shape(1);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
   // j-normal faces
      for( size_t j=0; j<nj; ++j )
     {
         for( size_t i=0; i<ni-1; ++i )
        {
            const Idx2<GT> ijl{i  ,j};
            const Idx2<GT> ijr{i+1,j};

            const auto acc_v = apply_stencil2( edge_func,
                                               ijl,ijr,
                                               src0, srcs... );

            dst(ijl) = accl_func( std::move(dst(ijl)),
                                  acc_v );

            dst(ijr) = accr_func( std::move(dst(ijr)),
                                  acc_v );
        }
     }

# ifdef _OPENMP
   # pragma omp parallel for
# endif
   // i-normal faces
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj-1; ++j )
        {
            const Idx2<GT> ijl{i,j  };
            const Idx2<GT> ijr{i,j+1};

            const auto acc_v = apply_stencil2( edge_func,
                                               ijl,ijr,
                                               src0, srcs... );

            dst(ijl) = accl_func( std::move(dst(ijl)),
                                  acc_v );

            dst(ijr) = accr_func( std::move(dst(ijr)),
                                  acc_v );
        }
     }

      return;
  }
}
