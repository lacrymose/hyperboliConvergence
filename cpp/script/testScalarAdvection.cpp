
# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <iostream>

   template<typename T>
   void f( T ){ std::cout << "f(typename)" << std::endl; }

   template<ImplementedVarSet VarT>
   void f( VarT ){ std::cout << "f(ImplementedVarSet)" << std::endl; }

   template<ScalarConservedVariables VarT> requires ImplementedVarSet<VarT>
   void f( VarT ){ std::cout << "f(ScalarConservedVariables)" << std::endl; }

   template<EulerConservedVariables VarT> requires ImplementedVarSet<VarT>
   void f( VarT ){ std::cout << "f(EulerConservedVariables)" << std::endl; }

   template<ScalarState StateT>
   void f( StateT ){ std::cout << "f(ScalarState)" << std::endl; };

   template<typename T>
   void g( T ){ std::cout << "g(typename)" << std::endl; }

   template<ImplementedVarSet VarT>
   void g( VarT ){ std::cout << "g(ImplementedVarSet)" << std::endl; }


   int main()
  {
//    constexpr LawType Law = LawType::ScalarAdvection;
      constexpr LawType Law = LawType::Euler;
      constexpr int nDim = 2;
      
      using BasisT = BasisType<Law>;

      using VarT = VariableSet<Law,nDim,BasisT::Conserved>;
      using StateT = State<Law,nDim>;

      using VarW = VariableSet<Law,nDim,BasisT::Characteristic>;

//    VarW w{4.,3.,2.};
      VarW w{4.,3.,2.,1.};
      VarT q;
      StateT state;
      StateT sl,sr;
      Species<Law> species;
      Geometry::Surface<nDim> face;

      std::cout << w[0] << std::endl;
      std::cout << w[1] << std::endl;
      std::cout << w[2] << std::endl;
      std::cout << std::endl;

      std::cout << q[0] << std::endl;
      std::cout << q[1] << std::endl;
      std::cout << q[2] << std::endl;
      std::cout << std::endl;

      CentralFlux<Law> central;

      FluxResult<Law,nDim> fr = central.flux( species, face, sl, sr );
      fr = central( species, face, sl, sr );
      fr = central( species, face,  q,  q );

      q+=fr.flux;

      std::cout << "true:   " << true  << std::endl;
      std::cout << "false:  " << false << std::endl;

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

      state = set2State( species, q );
      q = state2Set<VarT>( species, state );

      return 0;
  }
