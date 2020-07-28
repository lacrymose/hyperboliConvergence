
# pragma once

# include <utils/concepts.h>

# include <iostream>

   template<floating_point Real>
   struct UnsteadyTimeControls
  {
   // total number of timesteps
      size_t nTimesteps;

   // cfl number
      Real cfl;
  };

   template<floating_point Real>
   inline std::istream& operator>>( std::istream& is, UnsteadyTimeControls<Real>& tc )
  {
      is >> tc.nTimesteps >> tc.cfl;
      return is;
  }

