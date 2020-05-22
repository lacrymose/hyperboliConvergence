# ifndef FIELDOPS_H
# define FIELDOPS_H

# include <ode.h>
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

   template<typename SolutionType, typename ResidualType>
   void residual( const IdealGas2D::Species& gas,
                  const Array::Array1D< IdealGas2D::VariableSet<  SolutionType > >&   q0,
                  const Array::Array1D< IdealGas2D::VariableSet<  SolutionType > >&   q1,
                                        IdealGas2D::VariableDelta<ResidualType>&     res );

   template<typename SolutionType>
   void dualTimeResidual( const ODE::Implicit::MultiStep& otime,
                          const IdealGas2D::Species&        gas,
                          const Controls::GridControls1D&  grid,
                          const Types::Real                  dt,
                          const std::vector<Array::Array1D<IdealGas2D::VariableSet<SolutionType>>>&   q,
                          const std::vector<Array::Array1D<IdealGas2D::ConservedDelta>>&              r,
                                Array::Array1D<IdealGas2D::ConservedDelta>&                          rm );

   template<typename T>
   void centralExplicitSmoothing( const Controls::GridControls1D&  grid,
                                  const Types::Real               alpha,
                                  const Array::Array1D<T>&           r0,
                                        Array::Array1D<T>&           r1 );

   template<typename T>
   void centralImplicitSmoothing( const Controls::GridControls1D&  grid,
                                  const Types::Real               alpha,
                                  const Array::Array1D<T>&           r0,
                                        Array::Array1D<T>&           r1 );

   template<typename T>
   void multigridRestriction( const Controls::GridControls1D&  grid,
                              const Array::Array1D<T>&           rh,
                                    Array::Array1D<T>&          r2h );

   template<typename T>
   void multigridProlongation( const Controls::GridControls1D&  grid,
                               const Array::Array1D<T>&          e2h,
                                     Array::Array1D<T>&           eh );

# include <fieldOperations/gradientCalculation.ipp>
# include <fieldOperations/fluxResidual.ipp>
# include <fieldOperations/dualTimeResidual.ipp>

# include <fieldOperations/residual.ipp>

# include <fieldOperations/eulerForwardUpdate.ipp>

# include <fieldOperations/centralExplicitSmoothing.ipp>
# include <fieldOperations/centralImplicitSmoothing.ipp>

# include <fieldOperations/multigrid.ipp>

# endif
