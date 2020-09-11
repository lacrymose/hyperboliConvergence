
# include <solutionField/solutionField.h>

# include <conservationLaws/scalarAdvection/scalarAdvection.h>

# include <parallalg/array.h>

# include <vector>
# include <array>

# include <iostream>

constexpr LawType law = LawType::ScalarAdvection;
constexpr int nDim = 1;

using BasisT = BasisType<law>;
using Real = double;

constexpr BasisT basis = BasisT::Conserved;

   struct NonCopyable
  {
      NonCopyable() = default;
      NonCopyable( const NonCopyable&  ) = delete;
      NonCopyable(       NonCopyable&& ) = default;

      NonCopyable( int i ){};

      NonCopyable& operator=( const NonCopyable&  ) = delete;
      NonCopyable& operator=(       NonCopyable&& ) = default;

     ~NonCopyable() = default;
  };

   int main()
  {
//    std::array<NonCopyable,2> arr{0,1};
//    std::vector<NonCopyable>  vec(2,2);

      par::Shape<nDim> shape{4};

      SolutionField<law,nDim,basis,Real> q(shape);

      std::cout << q.q.shape(0)     << "\n";
      std::cout << q.qb.size()      << "\n";
      std::cout << q.qb[0].shape(0) << "\n";

      return 0;
  }
