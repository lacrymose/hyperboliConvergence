
# pragma once

# include <cassert>

// ---------- Species struct ----------

   auto get_species( ScalarAdvectionBases Basis )
  {
      assert( Basis == ScalarAdvectionBases::Conserved );

      return Species<LawType::ScalarAdvection,Real>{};
  }

   auto get_species( EulerBases Basis )
  {
      return get_air_species<Real>();
  }


// ---------- Scalar Advection initial states ----------

   auto initialLeft( ScalarAdvectionBases Basis )
  {
      using VarSet = VariableSet<LawType::ScalarAdvection,
                                 nDim,
                                 ScalarAdvectionBases::Conserved,
                                 Real>;
      constexpr Real u=1.;
      constexpr Real q=2.;
      return VarSet{u,q};
  }

   auto initialRight( ScalarAdvectionBases Basis )
  {
      using VarSet = VariableSet<LawType::ScalarAdvection,
                                 nDim,
                                 ScalarAdvectionBases::Conserved,
                                 Real>;
      constexpr Real u=1.;
      constexpr Real q=1.;
      return VarSet{u,q};
  }


// ---------- Euler initial states ----------

   auto initialLeft( EulerBases Basis )
  {
      using VarSet = VariableSet<LawType::Euler,
                                 nDim,
                                 EulerBases::Conserved,
                                 Real>;
      constexpr Real r =1.0;
      constexpr Real ru=0.0;
      constexpr Real re=2.5;
      return VarSet{ru,r,re};
  }
   
   auto initialRight( EulerBases Basis )
  {
      using VarSet = VariableSet<LawType::Euler,
                                 nDim,
                                 EulerBases::Conserved,
                                 Real>;
      constexpr Real r =0.125;
      constexpr Real ru=0.0;
      constexpr Real re=0.25;
      return VarSet{ru,r,re};
  }


// ---------- Write state to stream ----------

   void writeState( std::ofstream& os, const State<LawType::ScalarAdvection,1,Real>& s )
  {
      os << s.velocity(0) << " "
         << s.scalar()    << std::endl;
      return;
  }
   
   void writeState( std::ofstream& os, const State<LawType::Euler,1,Real>& s )
  {
      os << s.density()   << " "
         << s.velocity(0) << " "
         << s.pressure()  << " "
         << s.specificTotalEnthalpy()  << std::endl;
      return;
  }
 
