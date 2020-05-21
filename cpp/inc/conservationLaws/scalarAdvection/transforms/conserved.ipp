
   template<ScalarConservedVariables ConsVarT, ScalarState StateT, floating_point Real>
      requires   SameDim<ConsVarT,StateT>
              && SameFPType<ConsVarT,StateT,Real>
   ConsVarT state2Set( const Species<LawType::ScalarAdvection,Real>& species, const StateT& state )
  {
      constexpr int nd = dim_of_v<StateT>;

      ConsVarT qc;
      for( int i=0; i<nd; i++ )
     {
         qc[i]=state.velocity(i);
     }
      qc[nd]=state.scalar();

      return qc;
  }

   template<ScalarConservedVariables ConsVarT, floating_point Real>
      requires SameFPType<ConsVarT,Real>
   state_t<ConsVarT> set2State( const Species<LawType::ScalarAdvection,Real>& species, const ConsVarT& qc )
  {
      constexpr int nd = dim_of_v<ConsVarT>;
      using StateT = state_t<ConsVarT>;

      StateT state;
      for( int i=0; i<nd; i++ )
     {
         state.velocity(i)=qc[i];
     }
      state.scalar()=qc[nd];

      return state;
  }

