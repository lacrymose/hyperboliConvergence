
# pragma once

# include <conservationLaws/base/base.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

   template<typename ElemT>
   auto makeBoundaryArray( const par::Shape<1>& s )
  {
      using ArrayT = par::Array<ElemT,1>;
      using BoundaryArrays = std::array<ArrayT,2>;

      par::Shape<1> s01{2};

      return BoundaryArrays{ArrayT(s01),ArrayT(s01)};
  }

   template<typename ElemT>
   auto makeBoundaryArray( const par::Shape<2>& s )
  {
      using ArrayT = par::Array<ElemT,2>;
      using BoundaryArrays = std::array<ArrayT,4>;

      par::Shape<2> s01{s[1],2};
      par::Shape<2> s23{s[0],2};

      return BoundaryArrays{ArrayT(s01),ArrayT(s01),
                            ArrayT(s23),ArrayT(s23)};
  }

   template<typename ElemT>
   auto makeBoundaryArray( const par::Shape<3>& s )
  {
      using ArrayT = par::Array<ElemT,3>;
      using BoundaryArrays = std::array<ArrayT,6>;

      par::Shape<3> s01{s[1],s[2],2};
      par::Shape<3> s23{s[0],s[2],2};
      par::Shape<3> s45{s[0],s[1],2};

      return BoundaryArrays{ArrayT(s01),ArrayT(s01),
                            ArrayT(s23),ArrayT(s23),
                            ArrayT(s45),ArrayT(s45)};
  }

/*
 * data structure holding the solution field for a conservation law.
 *    solution over the internal domain
 *    solution over the boundary
 *    boundary condition types
 *    boundary reference solutions
 */
   template<ImplementedVarSet VarSet,
            int                 nDim>
   struct SolutionField
  {
      const static unsigned int nBoundaries=2*nDim;
      const static LawType Law = law_of_v<VarSet>;

      using VarDel = vardelta_t<VarSet>;

   // array type for field
      using VarField = par::Array<VarSet,nDim>;

   // enum for specifying boundary condition
      using BCType = BoundaryType<Law>;

   // internal and boundary solution fields
      VarField                         interior;
      std::array<VarField,nBoundaries> boundary;

   // boundary condition types
      std::array<BCType,nBoundaries> bcTypes;

   // par::Array only supports move construction, so same must be for SolutionField
      SolutionField() = delete;
      SolutionField( const SolutionField&  ) = delete;
      SolutionField(       SolutionField&& ) = default;

   // par::Array only supports move assignment, so same must be for SolutionField
      SolutionField& operator=( const SolutionField&  ) = delete;
      SolutionField& operator=(       SolutionField&& ) = default;

   // must be initialised with shape of domain
      SolutionField( const par::Shape<nDim>& s )
                     : interior(s),
                       boundary(makeBoundaryArray<VarSet>(s)){}
  };


   template<ImplementedVarSet VarSet,
            int                 nDim>
   void copy(       SolutionField<VarSet,nDim>& dst,
              const SolutionField<VarSet,nDim>& src )
  {
      assert( dst.interior.shape() == src.interior.shape() );

      par::copy( dst.interior, src.interior );
      for( unsigned int i=0; i<src.nBoundaries; i++ )
     {
         par::copy( dst.boundary[i], src.boundary[i] );
         dst.bcTypes[i] = src.bcTypes[i];
     }
  }

   template<ImplementedVarSet VarSet,
            int                 nDim>
   SolutionField<VarSet,nDim> copy( const SolutionField<VarSet,nDim>& src )
  {
      SolutionField<VarSet,nDim> dst(src.interior.shape());
      copy( dst,src );
      return dst;
  }



