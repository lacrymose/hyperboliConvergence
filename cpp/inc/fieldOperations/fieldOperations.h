# ifndef FIELDOPS_H
# define FIELDOPS_H

# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>

   template<char C>
   void fluxResidual( const IdealGas2D::Species& gas, const char boundaryCondition,
                      const Array1D< IdealGas2D::VariableSet<C> >&        q,
                            Array1D< IdealGas2D::ConservedVariables >&    r,
                            float& lmax );

   template<char C>
   void eulerForwardUpdate( const IdealGas2D::Species& gas, float dt,
                            const Array1D< IdealGas2D::VariableSet<C> >&        q0,
                                  Array1D< IdealGas2D::VariableSet<C> >&        q,
                            const Array1D< IdealGas2D::ConservedVariables >&    r  );

# include <fieldOperations/fluxResidual.ipp>
# include <fieldOperations/eulerForwardUpdate.ipp>

# endif
