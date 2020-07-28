
# pragma once

# include <ode.h>

# include <conservationLaws/base/base.h>

# include <parallalg/array.h>

# include <vector>

   template<LawType Law, floating_point Real>
   par::Array<FluxResult<Law,1,Real>,1>
      rungeKuttaAccumulation( const ODE::Explicit::RungeKutta<Real>& rk,
                              const unsigned int        stg,
                              const std::vector<par::Array<FluxResult<Law,1,Real>,1>>& resStage )
  {
      assert( resStage.size()<stg );

      using FluxRes = FluxResult<Law,1,Real>;

      const par::Shape<1> shape(resStage[0].shape());
      par::Array<FluxRes,1> resTotal(shape);

      par::fill( resTotal, FluxRes{} );

      const size_t ni=shape[0];

      for( size_t i=0; i<ni; i++ )
     {
         for( unsigned int l=0; l<=stg; l++ )
        {
            resTotal[{i}].flux += rk.alpha[stg][l]*resStage[l][{i}].flux;
        }
     }
      return resTotal;
  }

   template<LawType Law, floating_point Real>
   par::Array<FluxResult<Law,2,Real>,2>
      rungeKuttaAccumulation( const ODE::Explicit::RungeKutta<Real>& rk,
                              const unsigned int        stg,
                              const std::vector<par::Array<FluxResult<Law,2,Real>,2>>& resStage )
  {
      assert( resStage.size()<stg );

      using FluxRes = FluxResult<Law,2,Real>;

      const par::Shape<2> shape(resStage[0].shape());
      par::Array<FluxRes,2> resTotal(shape);

      par::fill( resTotal, FluxRes{} );

      const size_t ni=shape[0];
      const size_t nj=shape[1];

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
      return resTotal;
  }

   template<LawType Law, floating_point Real>
   par::Array<FluxResult<Law,3,Real>,3>
      rungeKuttaAccumulation( const ODE::Explicit::RungeKutta<Real>& rk,
                              const unsigned int        stg,
                              const std::vector<par::Array<FluxResult<Law,3,Real>,3>>& resStage )
  {
      assert( resStage.size()<stg );

      using FluxRes = FluxResult<Law,3,Real>;

      const par::Shape<3> shape(resStage[0].shape());
      par::Array<FluxRes,3> resTotal(shape);

      par::fill( resTotal, FluxRes{} );

      const size_t ni=shape[0];
      const size_t nj=shape[1];
      const size_t nk=shape[2];

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
      return resTotal;
  }

