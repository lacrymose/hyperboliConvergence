
# pragma once

# include <conservationLaws/base/base.h>

# include <parallalg/array.h>

   template<typename ElemT>
   auto makeBoundaryArray( const par::Shape<1>& s )
  {
      using ArrayT = par::Array<ElemT,1>;
      using BoundaryArrays = std::array<ArrayT,2>;

      par::Shape<1> s01({2});

      return BoundaryArrays{ArrayT(s01),ArrayT(s01)};
  }

   template<typename ElemT>
   auto makeBoundaryArray( const par::Shape<2>& s )
  {
      using ArrayT = par::Array<ElemT,2>;
      using BoundaryArrays = std::array<ArrayT,4>;

      par::Shape<2> s01({s[1],2});
      par::Shape<2> s23({s[0],2});

      return BoundaryArrays{ArrayT(s01),ArrayT(s01),
                            ArrayT(s23),ArrayT(s23)};
  }

   template<typename ElemT>
   auto makeBoundaryArray( const par::Shape<3>& s )
  {
      using ArrayT = par::Array<ElemT,3>;
      using BoundaryArrays = std::array<ArrayT,6>;

      par::Shape<3> s01({s[1],s[2],2});
      par::Shape<3> s23({s[0],s[2],2});
      par::Shape<3> s45({s[0],s[1],2});

      return BoundaryArrays{ArrayT(s01),ArrayT(s01),
                            ArrayT(s23),ArrayT(s23),
                            ArrayT(s45),ArrayT(s45),};
  }

/*
 * data structure holding the solution field for a conservation law.
 *    solution over the internal domain
 *    solution over the boundary
 *    boundary condition types
 *    boundary reference solutions
 */
   template<LawType          Law,
            int             nDim,
            BasisType<Law> Basis,
            floating_point  Real>
   struct SolutionField
  {
      const static unsigned int nBoundaries=2*nDim;

      using VarSet = VariableSet<  Law,nDim,Basis,Real>;
      using VarDel = VariableDelta<Law,nDim,Basis,Real>;

      using VarField = par::Array<VarSet,nDim>;

      using BCType = BoundaryType<Law>;

      par::Shape<nDim> shape;

      VarField   q;
      std::array<VarField,nBoundaries>     qb;
      std::array<BCType,  nBoundaries> bcType;

   // par::Array only supports move construction, so same must be for SolutionField
      SolutionField() = delete;
      SolutionField( const SolutionField&  ) = delete;
      SolutionField(       SolutionField&& ) = default;

   // par::Array only supports move assignment, so same must be for SolutionField
      SolutionField& operator=( const SolutionField&  ) = delete;
      SolutionField& operator=(       SolutionField&& ) = default;

      SolutionField( const par::Shape<nDim>& s )
                     : shape(s),
                       q(s),
                       qb(makeBoundaryArray<VarSet>(s)){}
  };

