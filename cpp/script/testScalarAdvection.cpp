
# include <conservationLaws/scalarAdvection/scalarAdvection.h>

# include <geometry/geometry.h>

# include <iostream>

   template<typename T>
   void f( T ){ std::cout << "f(typename)" << std::endl; }

   template<ImplementedVarSet VarT>
   void f( VarT ){ std::cout << "f(ImplementedVarSet)" << std::endl; }

   template<ScalarConservedVariables VarT> requires ImplementedVarSet<VarT>
   void f( VarT ){ std::cout << "f(ScalarConservedVariables)" << std::endl; }

   template<ScalarState StateT>
   void f( StateT ){ std::cout << "f(ScalarState)" << std::endl; }

   template<typename T>
   void g( T ){ std::cout << "g(typename)" << std::endl; }

   template<ImplementedVarSet VarT>
   void g( VarT ){ std::cout << "g(ImplementedVarSet)" << std::endl; }


   int main()
  {
      using FPType = float;
      constexpr LawType Law = LawType::ScalarAdvection;
      constexpr int nDim = 2;
      
      using BasisT = BasisType<Law>;

      using VarT = VariableSet<Law,nDim,BasisT::Conserved,FPType>;
      using StateT = State<Law,nDim,FPType>;
      using SpeciesT = Species<Law,FPType>;

      VarT q{FPType(4.),FPType(3.),FPType(2.)};
      StateT state{};
//    StateT sl{},sr{};
      SpeciesT species{};
//    Geometry::Surface<nDim,FPType> face{};


      std::cout << "true:   " << true  << std::endl;
      std::cout << "false:  " << false << std::endl;
      std::cout << std::endl;

      std::cout << std::is_same_v< state_t<VarT>,
                                   StateT> << std::endl;

      std::cout << std::is_same_v< species_t<VarT>,
                                   SpeciesT> << std::endl;
      std::cout << std::endl;


      std::cout << ScalarState<StateT> << std::endl;
      std::cout << ScalarConservedVariables<VarT> << std::endl;
      std::cout << std::endl;

      std::cout << std::is_same_v<fptype_of_t<StateT>,FPType> << std::endl;
      std::cout << std::is_same_v<fptype_of_t<VarT>,FPType> << std::endl;
      std::cout << SameFPType<StateT,VarT,FPType> << std::endl;
      std::cout << std::endl;

//    state = set2State( species, q );
//    q = state2Set<VarT>( species, state );

//    CentralFlux<Law> central;

/*
      FluxResult<Law,nDim> fr = central.flux( species, face, sl, sr );
      fr = central( species, face, sl, sr );
      fr = central( species, face,  q,  q );

      q+=fr.flux;

      std::cout << is_State_v<StateT> << std::endl;
      std::cout << has_same_law_v<StateT,law_constant<Law>> << std::endl;
      std::cout << ScalarState<StateT> << std::endl;

      std::cout << is_specialised_VarSet_v<VarT,Law,BasisT::Conserved> << std::endl;
      std::cout << std::endl;

      std::cout << has_same_dim_v<VarT,StateT> << std::endl;
      std::cout << SameDim<VarT,StateT> << std::endl;
      std::cout << std::endl;

      std::cout << FluxFunctor<CentralFlux<Law>,Law> << std::endl;
      std::cout << std::endl;

      f(w);
      f(q);
      f(state);

      std::cout << std::endl;

      g(w);
      g(q);


      for( int i=0; i<nDim+1; i++ ){ q[i]=i; }
*/

      return 0;
  }
