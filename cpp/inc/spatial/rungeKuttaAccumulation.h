
# pragma once

# include <ode.h>

# include <conservationLaws/base/base.h>

# include <parallalg/array.h>

# include <vector>


   template<LawType Law, int nDim, floating_point Real>
   par::Array<FluxResult<Law,nDim,Real>,nDim>
      rungeKuttaAccumulation( const ODE::Explicit::RungeKutta<Real>& rk,
                              const unsigned int        stg,
                              const std::vector<par::Array<FluxResult<Law,nDim,Real>,nDim>>& resStage )
  {
      using FluxRes = FluxResult<Law,nDim,Real>;
      par::Array<FluxRes,nDim> resTotal(resStage[0].shape());
      rungeKuttaAccumulation( rk, stg, resStage, resTotal );
      return resTotal;
  }

   template<LawType Law, floating_point Real>
   void rungeKuttaAccumulation( const ODE::Explicit::RungeKutta<Real>& rk,
                                const unsigned int        stg,
                                const std::vector<par::Array<FluxResult<Law,1,Real>,1>>& resStage,
                                      par::Array<FluxResult<Law,1,Real>,1>&              resTotal )
  {
      assert( resStage.size()==rk.nstages );
      for( const auto& r : resStage ){ assert( r.shape()==resTotal.shape() ); }

      using FluxRes = FluxResult<Law,1,Real>;
      par::fill( resTotal, FluxRes{} );

      const size_t ni=resTotal.shape(0);

      for( size_t i=0; i<ni; i++ )
     {
         for( unsigned int l=0; l<=stg; l++ )
        {
            resTotal[{i}].flux += rk.alpha[stg][l]*resStage[l][{i}].flux;
        }
     }
  }

   template<LawType Law, floating_point Real>
   void rungeKuttaAccumulation( const ODE::Explicit::RungeKutta<Real>& rk,
                                const unsigned int        stg,
                                const std::vector<par::Array<FluxResult<Law,2,Real>,2>>& resStage,
                                      par::Array<FluxResult<Law,2,Real>,2>&              resTotal )
  {
      assert( resStage.size()==rk.nstages );
      for( const auto& r : resStage ){ assert( r.shape()==resTotal.shape() ); }

      using FluxRes = FluxResult<Law,2,Real>;
      par::fill( resTotal, FluxRes{} );

      const size_t ni=resTotal.shape(0);
      const size_t nj=resTotal.shape(1);

      for( size_t i=0; i<ni; i++ )
     {
         for( size_t j=0; j<nj; j++ )
        {
            for( unsigned int l=0; l<=stg; l++ )
           {
               resTotal[{i,j}].flux += rk.alpha[stg][l]*resStage[l][{i,j}].flux;
           }
        }
     }
  }

   template<LawType Law, floating_point Real>
   void rungeKuttaAccumulation( const ODE::Explicit::RungeKutta<Real>& rk,
                                const unsigned int        stg,
                                const std::vector<par::Array<FluxResult<Law,3,Real>,3>>& resStage,
                                      par::Array<FluxResult<Law,3,Real>,3>&              resTotal )
  {
      assert( resStage.size()==rk.nstages );
      for( const auto& r : resStage ){ assert( r.shape()==resTotal.shape() ); }

      using FluxRes = FluxResult<Law,3,Real>;
      par::fill( resTotal, FluxRes{} );

      const size_t ni=resTotal.shape(0);
      const size_t nj=resTotal.shape(1);
      const size_t nk=resTotal.shape(2);

      for( size_t i=0; i<ni; i++ )
     {
         for( size_t j=0; j<nj; j++ )
        {
            for( size_t k=0; k<nk; k++ )
           {
               for( unsigned int l=0; l<=stg; l++ )
              {
                  resTotal[{i,j,k}].flux += rk.alpha[stg][l]*resStage[l][{i,j,k}].flux;
              }
           }
        }
     }
  }

