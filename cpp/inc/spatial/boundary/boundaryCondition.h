
# pragma once

# include <conservationLaws/base/base.h>

# include <utils/concepts.h>
# include <utils/utils.h>

# include <array>
# include <tuple>
# include <type_traits>

   template<LawType              Law,
            BoundaryType<Law> BCType,
            typename...        Funcs>
   struct BoundaryCondition
  {
      std::tuple<Funcs...> funcs;

      BoundaryCondition( const Funcs&... fs ) : funcs(fs...){};

      constexpr BoundaryType<Law> type() const { return BCType; }

      template<size_t N>
         requires (N < sizeof...(Funcs))
      const auto& func() const { return std::get<N>(funcs); }
  };

   template<typename              BCType,
            typename...    BoundaryConds,
            typename    BoundaryFunction,
            typename...             Args>
   void selectBoundaryCondition( const BCType               typeToMatch,
                                 const std::tuple<BoundaryConds...> bcs,
                                 const BoundaryFunction&        bc_func,
                                       Args&&...                   args )
  {
   // compare on boundary condition type
      const auto get_type = []( const auto& bc ){ return bc.type(); };
   // check if equal
      const auto is_equal = []( auto l, auto r ){ return l==r; };

      utils::selectFromTuple<0>( typeToMatch, bcs,
                                 get_type, is_equal, bc_func,
                                 std::forward<Args>(args)... );
      return;
  }

   template<LawType              Law,
            BoundaryType<Law> BCType,
            typename...        Funcs>
   auto make_BCond( Funcs&&... funcs )
  {
      return BoundaryCondition<Law,BCType,Funcs...>{funcs...};
  }

   template<LawType Law>
   auto make_periodic_BCond()
  {
      return make_BCond<Law,BoundaryType<Law>::Periodic>();
  }

   template<LawType              Law,
            BoundaryType<Law> BCType>
   auto make_boundary_update_lambda()
  {
      return [] <int                    nDim,
                 floating_point         Real,
                 ImplementedVarSet   SolVarT>
               ( const Species<Law,Real>&            spc,
                 const geom::Surface<nDim,Real>&    face,
                 const geom::Volume<nDim,Real>&    celli,
                 const SolVarT&                qinterior,
                 const SolVarT&                qboundary ) -> SolVarT
//       requires ConsistentTypes<Law,nDim,Real,SolVarT>
     {
         using tag = std::integral_constant<BoundaryType<Law>,BCType>;
         return boundaryUpdate( tag{}, spc, face, celli, qinterior, qboundary );
     };
  }

   template<LawType              Law,
            BoundaryType<Law> BCType>
   auto make_boundary_flux_lambda()
  {
      return [] <int                    nDim,
                 floating_point         Real,
                 ImplementedVarSet   SolVarT,
                 ImplementedVarSet   SolDelT>
               ( const Species<Law,Real>&            spc,
                 const geom::Surface<nDim,Real>&    face,
                 const geom::Volume<nDim,Real>&    celli,
                 const SolDelT                        dq,
                 const SolVarT&                qinterior,
                 const SolVarT&                qboundary,
                 const SolVarT&                qfarfield ) -> FluxResult<Law,nDim,Real>
//       requires ConsistentTypes<Law,nDim,Real,SolVarSet,SolDelT>
     {
         using tag = std::integral_constant<BoundaryType<Law>,BCType>;
         return boundaryFlux( tag{}, spc, face, celli, dq, qinterior, qboundary, qfarfield );
     };
  }

   template<LawType              Law,
            BoundaryType<Law> BCType>
   auto make_ghostCell_BCond()
  {
      return make_BCond<Law,BCType>(make_boundary_update_lambda<Law,BCType>());
  }

   template<LawType              Law,
            BoundaryType<Law> BCType>
   auto make_flux_BCond()
  {
      return make_BCond<Law,BCType>(make_boundary_update_lambda<Law,BCType>(),
                                    make_boundary_flux_lambda<Law,BCType>());
  }

/*
   template<typename  T,
            LawType Law,
            ...>
   concept bool BoundaryConditionUpdate = 
      std::is_invocable_r<

   template<typename T,
            ...>
   concept bool BoundaryConditionFlux = 
*/
