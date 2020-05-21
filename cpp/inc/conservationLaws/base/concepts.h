
# pragma once

# include <conservationLaws/base/type-traits.h>
# include <conservationLaws/base/declarations.h>

# include <geometry/geometry.h>

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
 * use the same floating point precision
 */
   template<typename... Ts>
   concept bool SameFPType = has_same_fptype_v<Ts...>;

/*
 * dimension is physically realistic (1, 2 or 3)
 */
   template<typename T>
   concept bool HasValidDim = has_valid_dim_v<T>;


/*
 * The minimal functionality has been implemented for this VariableSet<LawType,int,Basis>
 *    Transformations to/from a State is defined
 */
   template<typename T>
   concept bool ImplementedVarSet = 
      is_VariableSet_v<T>
   && HasValidDim<T>
   && requires( T q, species_t<T> species, state_t<T> state )
     {
         { set2State(    species,     q ) } -> state_t<T>;
         { state2Set<T>( species, state ) } -> T;
     };

/*
 * The minimal functionality has been implemented for this VariableSet<LawType,int,Basis>
 *    Transformations to/from a State is defined
 */
   template<typename T>
   concept bool ImplementedVarDelta = 
      is_VariableDelta_v<T>
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
      ImplementedVarSet<  VariableSet<  Law,2,BasisType<Law>::Conserved,double>>
   && ImplementedVarDelta<VariableDelta<Law,2,BasisType<Law>::Conserved,double>>
   && requires( Species<Law,double> species, State<Law,2,double> state, Geometry::Direction<2,double> normal, int ns, int nv )
     {
         ns = nScalarQuantities<Law>;
         nv = nVectorQuantities<Law>;

         species;
         state;

         { exactFlux( species, normal, state ) } -> FluxResult<Law,2,double>;
     };

/*
 * A FluxFunctor is useable as an interface flux between cells in the Finite Volume Method
 *    simplest way of satisfying this concept is to combine a type satisfying the FluxImplementation concept with the CRTP base type FluxInterface
 */
   template<typename T, LawType Law>
   concept bool FluxFunctor =
      ImplementedLawType<Law>
   && requires( T Flux, Species<Law,double> species,  Geometry::Surface<2,double> face,
                State<Law,2,double> sl, State<Law,2,double> sr,
                VariableSet<Law,2,BasisType<Law>::Conserved,double> ql,
                VariableSet<Law,2,BasisType<Law>::Conserved,double> qr )
     {
         { Flux( species, face, sl, sr ) } ->  FluxResult<Law,2,double>;
         { Flux( species, face, ql, qr ) } ->  FluxResult<Law,2,double>;
     };


