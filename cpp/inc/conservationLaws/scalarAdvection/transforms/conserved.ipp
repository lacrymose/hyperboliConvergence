
   template<ScalarConservedVariables ConsVarT, ScalarState StateT>
      requires SameDim<ConsVarT,StateT>
   ConsVarT state2Set( const Species<LawType::ScalarAdvection>& species, const StateT& state )
  {
      constexpr int nd = dim_of_v<StateT>;

      ConsVarT qc;
      for( int i=0; i<nd; i++ )
     {
         qc[i]=state.velocity(i);
     }
      qc[nd]=species.scale*state.scalar();

      return qc;
  }

   template<ScalarConservedVariables ConsVarT>
   state_t<ConsVarT> set2State( const Species<LawType::ScalarAdvection>& species, const ConsVarT& qc )
  {
      constexpr int nd = dim_of_v<ConsVarT>;
      using StateT = state_t<ConsVarT>;

      StateT state;
      for( int i=0; i<nd; i++ )
     {
         state.velocity(i)=qc[i];
     }
      state.scalar()=qc[nd]/species.scale;

      return state;
  }

