
# pragma once

# include <conservationLaws/base/type-traits.h>
# include <conservationLaws/base/declarations.h>

# include <utils/concepts.h>

// concepts

/*
 * all types exist in the same spatial dimensions
 */
   template<typename... Ts>
   concept bool SameDim = has_same_dim_v<Ts...>;

/*
 * all types relate to the same conservation law
 */
   template<typename... Ts>
   concept bool SameLaw = has_same_law_v<Ts...>;

/*
 * dimension is physically realistic (1, 2 or 3)
 */
   template<typename T>
   concept HasValidDim = has_valid_dim_v<T>;

/*
 * The minimal functionality has been implemented for this VariableSet<LawType,int,Basis>
 *    Transformations to/from a State is defined
 */
   template<typename T>
   concept bool ImplementedVarSet = 
      same_template_as<T,VariableSet>
   && HasValidDim<T>
   && requires( T q, species_t<T> species, state_t<T> state )
     {
         { transform<state_t<T>>( species,     q ) } -> state_t<T>;
         { transform<T>(          species, state ) } -> T;
     };

/*
 * The minimal functionality has been implemented for this VariableSet<LawType,int,Basis>
 *    Transformations to/from a State is defined
 */
   template<typename T>
   concept bool ImplementedVarDelta = 
      same_template_as<T,VariableDelta>
   && HasValidDim<T>;

/*
 * The minimal functionality has been implemented for this LawType
 *    The dimension of the phase space is defined
 *    Species and State specialisations are defined
 *    The exact physical flux function is defined
 *    Conserved variables are minimally implemented
 */
   template<LawType Law>
   concept bool ImplementedLawType = 
      ImplementedVarSet<  VariableSet<  Law,2,BasisType<Law>::Conserved>>
   && ImplementedVarDelta<VariableDelta<Law,2,BasisType<Law>::Conserved>>
   && requires( Species<Law> species, State<Law,2> state, Geometry::Surface<2> face )
     {
         ns = nScalarQuantities<Law>;
         nv = nVectorQuantities<Law>;

         species;
         state;

         { exactFlux( species, face, state ) } -> FluxResult<Law,2>;
     };

/*
 * A FluxImplementation provides the minimum functionality to be used as an interface flux.
 *    Generally a type satisfying FluxImplementation will be used with the CRTP class FluxInterface, together satisfying the FluxFunctor concept
 */
   template<typename T, LawType Law>
   concept bool FluxImplementation =
      ImplementedLawType<Law>
   && requires( T Flux, Species<Law> species,  Geometry::Surface<2> face,
                State<Law,2> sl, State<Law,2> sr )
     {
         { Flux.flux( species, face, sl, sr ) } ->  FluxResult<Law,2>;
     };

/*
 * A FluxFunctor is useable as an interface flux between cells in the Finite Volume Method
 *    simplest way of satisfying this concept is to combine a type satisfying the FluxImplementation concept with the CRTP base type FluxInterface
 */
   template<typename T, LawType Law>
   concept bool FluxFunctor =
      ImplementedLawType<Law>
   && requires( T Flux, Species<Law> species,  Geometry::Surface<2> face,
                State<Law,2> sl, State<Law,2> sr,
                VariableSet<Law,2,BasisType<Law>::Conserved> ql,
                VariableSet<Law,2,BasisType<Law>::Conserved> qr )
     {
         { Flux( species, face, sl, sr ) } ->  FluxResult<Law,2>;
         { Flux( species, face, ql, qr ) } ->  FluxResult<Law,2>;
     };


