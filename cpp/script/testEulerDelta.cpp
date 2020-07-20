
# include <conservationLaws/euler/euler.h>

# include <limits>

# include <iostream>
# include <iomanip>

using Real = long double;

constexpr LawType Law = LawType::Euler;

constexpr int nDim = 2;

using BasisT = BasisType<Law>;

constexpr BasisT ConsBasis = BasisT::Conserved;
constexpr BasisT PrimBasis = BasisT::Primitive;
constexpr BasisT ViscBasis = BasisT::Viscous;
constexpr BasisT CharBasis = BasisT::Characteristic;

using ConsVarT = VariableSet<Law,nDim,ConsBasis,Real>;
using PrimVarT = VariableSet<Law,nDim,PrimBasis,Real>;
using ViscVarT = VariableSet<Law,nDim,ViscBasis,Real>;

using ConsDelT = VariableDelta<Law,nDim,ConsBasis,Real>;
using PrimDelT = VariableDelta<Law,nDim,PrimBasis,Real>;
using ViscDelT = VariableDelta<Law,nDim,ViscBasis,Real>;
using CharDelT = VariableDelta<Law,nDim,CharBasis,Real>;

using StateT = State<Law,nDim,Real>;

   int main()
  {
      std::cout << std::scientific;
      std::cout.precision(8);

      const Species<Law,Real> gas = get_air_species<Real>();

      const ViscVarT qv0{Real(30.),Real(50.),Real(273.),Real(1.e5)};
      const PrimVarT qp0 = set2Set<PrimVarT>( gas, qv0 );

      const ViscDelT dqv0{Real(1.),Real(2.),Real(3.),  Real(100.)};
      const PrimDelT dqp0{Real(1.),Real(2.),Real(0.05),Real(100.)};

      const ViscDelT dqv1 = delta2Delta<ViscDelT>( gas, qv0, dqp0 );
      const PrimDelT dqp1 = delta2Delta<PrimDelT>( gas, qp0, dqv0 );

      const ViscDelT dqv2 = delta2Delta<ViscDelT>( gas, qp0, dqp1 );
      const PrimDelT dqp2 = delta2Delta<PrimDelT>( gas, qv0, dqv1 );

      std::cout << "dqv0\n";
      std::cout << dqv0 << std::endl;
      std::cout << "dqv2\n";
      std::cout << dqv2 << std::endl;
      std::cout << "diff\n";
      std::cout << dqv0-dqv2 << std::endl;

      std::cout << std::endl;

      std::cout << "dqp0\n";
      std::cout << dqp0 << std::endl;
      std::cout << "dqp2\n";
      std::cout << dqp2 << std::endl;
      std::cout << "diff\n";
      std::cout << dqp0-dqp2 << std::endl;

      std::cout << std::endl;

      std::cout << "dqp1\n";
      std::cout << dqp1 << std::endl;
      std::cout << "dqv1\n";
      std::cout << dqv1 << std::endl;

      std::cout << std::endl;

      std::cout << "characteristic:\n";

      const ViscDelT dqv3 = delta2Delta<ViscDelT>( gas, qv0,
                                                   delta2Delta<CharDelT>( gas, qv0, dqv0 ) );

      std::cout << "dqv0\n";
      std::cout << dqv0 << std::endl;
      std::cout << "dqv3\n";
      std::cout << dqv3 << std::endl;
      std::cout << "diff\n";
      std::cout << dqv0-dqv3 << std::endl;

      return 0;
  }
