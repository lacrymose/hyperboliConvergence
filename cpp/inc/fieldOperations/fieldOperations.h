# ifndef FIELDOPS_H
# define FIELDOPS_H

# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>

   template<typename SolutionType, typename Flux>
   void fluxResidual( const IdealGas2D::Species& gas, const char boundaryCondition, Flux flux,
                            const Array::Array1D< IdealGas2D::VariableSet< SolutionType > >& q,
                                  Array::Array1D< IdealGas2D::ConservedDelta >&              r,
                                                                         float&           lmax  );

   template<typename SolutionType>
   void eulerForwardUpdate( const IdealGas2D::Species& gas, float dt,
                            const Array::Array1D< IdealGas2D::VariableSet< SolutionType > >& q0,
                                  Array::Array1D< IdealGas2D::VariableSet< SolutionType > >&  q,
                            const Array::Array1D< IdealGas2D::ConservedDelta >&               r,
                                                    IdealGas2D::ConservedDelta&             res  );

# include <fieldOperations/fluxResidual.ipp>
# include <fieldOperations/eulerForwardUpdate.ipp>

# endif
