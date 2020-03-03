# ifndef FIELDOPS_H
# define FIELDOPS_H

# include <controls.h>
# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>

# include <types.h>

   template<typename SolutionType>
   void gradientCalculation( const Controls::GridControls1D&   grid,
                             const Array::Array1D< IdealGas2D::VariableSet<   SolutionType > >&  q,
                                   Array::Array1D< IdealGas2D::VariableDelta< SolutionType > >& dq );

   template<typename SolutionType, typename Flux, typename Limiter>
   void fluxResidual( const IdealGas2D::Species& gas, const char boundaryCondition,
                      const Flux flux, const Limiter limiter,
                      const Array::Array1D< IdealGas2D::VariableSet<   SolutionType > >&  q,
                      const Array::Array1D< IdealGas2D::VariableDelta< SolutionType > >& dq,
                            Array::Array1D< IdealGas2D::ConservedDelta >&                 r,
                                                             Types::Real&              lmax  );

   template<typename SolutionType>
   void eulerForwardLinearUpdate( const IdealGas2D::Species& gas, Types::Real dt,
                                  const Array::Array1D< IdealGas2D::VariableSet< SolutionType > >& q0,
                                        Array::Array1D< IdealGas2D::VariableSet< SolutionType > >&  q,
                                  const Array::Array1D< IdealGas2D::ConservedDelta >&               r,
                                                        IdealGas2D::ConservedDelta&               res  );

   template<typename SolutionType>
   void eulerForwardNonlinearUpdate( const IdealGas2D::Species& gas, Types::Real dt,
                                     const Array::Array1D< IdealGas2D::VariableSet< SolutionType > >& q0,
                                           Array::Array1D< IdealGas2D::VariableSet< SolutionType > >&  q,
                                     const Array::Array1D< IdealGas2D::ConservedDelta >&               r,
                                                           IdealGas2D::ConservedDelta&               res  );

# include <fieldOperations/gradientCalculation.ipp>
# include <fieldOperations/fluxResidual.ipp>
# include <fieldOperations/eulerForwardUpdate.ipp>

# endif
