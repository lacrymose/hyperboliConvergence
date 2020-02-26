# ifndef TIMESTEPPING_H
# define TIMESTEPPING_H

# include <ode.h>
# include <controls.h>
# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>
# include <fieldOperations/fieldOperations.h>

# include <types.h>

# include <iostream>

namespace TimeStepping
{
   template<typename Flux, typename SolutionType>
   void ExplicitEuler( const Controls::TimeSteppingControls& time,
                       const Controls::GridControls1D&       grid,
                       const IdealGas2D::Species&             gas,
                       const Flux                            flux,
                             Array::Array1D< IdealGas2D::VariableSet< SolutionType > >& q,
                       const bool print );
}

# include <timestepping/explicitEuler.ipp>

# endif
