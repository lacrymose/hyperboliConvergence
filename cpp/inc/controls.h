# ifndef CONTROLS_H
# define CONTROLS_H

# include <iostream>

namespace Controls
{
   struct GridControls1D
  {
   // number of grid cells
      int n;

   // type of boundary condition
      char boundaryCondition;
  };

   inline std::istream& operator>>( std::istream& is, GridControls1D& g )
  {
      is >> g.n >> g.boundaryCondition;
      return is;
  }

   struct TimeSteppingControls
  {
   // total number of timesteps
      int nt;

   // cfl number
      float cfl;
  };

   inline std::istream& operator>>( std::istream& is, TimeSteppingControls& t )
  {
      is >> t.nt >> t.cfl;
      return is;
  }
}

# endif
