
# include <idealGas2D/idealGas2D.h>
# include <array1D/array1D.h>

# include <iostream>

   int main()
  {
      IdealGas2D::Species  gas;
      IdealGas2D::VariableSet<IdealGas2D::VariableType<'c'>>     qc0;
      IdealGas2D::VariableDelta<IdealGas2D::VariableType<'c'>>  dqc0,dqc1;

      gas.air();

      qc0=1.8;
      std::cout << qc0 << std::endl;

      dqc0=1.8;
      dqc0/=1.35;
      std::cout << dqc0 << std::endl;

      dqc0*=0.75;
      std::cout << dqc0 << std::endl;

      dqc1=3.2;
      qc0+=dqc1;
      std::cout << qc0 << std::endl;

      qc0-=dqc0;
      std::cout << qc0 << std::endl;

/*
      qc2 = qc0+qc1;
      std::cout << qc2 << std::endl;

      qc2 = qc0-qc1;
      std::cout << qc2 << std::endl;

      qc2 = qc0*2;
      std::cout << qc2 << std::endl;

      qc2 = 2*qc0;
      std::cout << qc2 << std::endl;

      qc2 = qc0/2;
      std::cout << qc2 << std::endl;

      IdealGas2D::VariableSet<'c'> qc3(gas,qc2);
      std::cout << qc3 << std::endl;

      IdealGas2D::VariableSet<'v'> qv0;
      std::cout << qv0 << std::endl;

      qv0=10.;
      IdealGas2D::VariableSet<'c'> qc4(gas,qv0);
      std::cout << qc4 << std::endl;

      int n=3;
      Array1D<IdealGas2D::ConservedVariables> arrayCV(n);
      Array1D<IdealGas2D::ConservedVariables> arrayCW(n);

      arrayCV = 6.7;
      for( int i=0; i<n; i++ )
     {
         std::cout << arrayCV[i] << std::endl;
     }

      arrayCV*= 0.5;
      for( int i=0; i<n; i++ )
     {
         std::cout << arrayCV[i] << std::endl;
     }

      arrayCW = 10.5;
      for( int i=0; i<n; i++ )
     {
         std::cout << arrayCW[i] << std::endl;
     }

      arrayCV+= arrayCW;
      for( int i=0; i<n; i++ )
     {
         std::cout << arrayCV[i] << std::endl;
     }

      arrayCW = arrayCW*0.1;
      for( int i=0; i<n; i++ )
     {
         std::cout << arrayCW[i] << std::endl;
     }

      arrayCW = arrayCW-arrayCV;
      for( int i=0; i<n; i++ )
     {
         std::cout << arrayCW[i] << std::endl;
     }

*/
      return 0;
  }
