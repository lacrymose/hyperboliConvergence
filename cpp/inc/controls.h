
# pragma once

# include <utils/concepts.h>

# include <iostream>

/*
 * hold parameters for controlling an time-accurate time marching routine
 */
   template<floating_point Real>
   struct UnsteadyTimeControls
  {
   // total number of timesteps
      size_t nTimesteps;

   // cfl number
      Real cfl;
  };

   template<floating_point Real>
   std::istream& operator>>( std::istream& is, UnsteadyTimeControls<Real>& tc )
  {
      is >> tc.nTimesteps >> tc.cfl;
      return is;
  }

/*
 * hold parameters for controlling an steady-state convergence time marching routine
 */
   template<floating_point Real>
   struct SteadyTimeControls
  {
   // total number of timesteps
      size_t nTimesteps;

   // cfl number
      Real cfl;

   // residual drop convergence criteria
      Real residual_drop;

   // residual smoothing?
      bool use_residual_smoothing;

   // residual smoothing cfl factor
      Real smoothing_factor;
  };

