
# include <idealGas2D/test-VariableDelta.h>
# include <types.h>

   void Test_IdealGas2D_VariableDelta::setUp()
  {
      dq0=0.;
      dq1=0.;
  }

   void Test_IdealGas2D_VariableDelta::tearDown(){}

   void Test_IdealGas2D_VariableDelta::test_constructors()
  {
   // default constructor
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0.var[0], 0, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0.var[1], 0, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0.var[2], 0, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0.var[3], 0, Types::EPS );

      dq0.var[0]=1.;
      dq0.var[1]=2.;
      dq0.var[2]=3.;
      dq0.var[3]=4.;

   // copy constructor
      dq1=IdealGas2D::VariableDelta<IdealGas2D::VariableType<'0'>>(dq0);

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1.var[0], dq0.var[0], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1.var[1], dq0.var[1], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1.var[2], dq0.var[2], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1.var[3], dq0.var[3], Types::EPS );

   // copy constructor as dummy for linear jacobian (which needs a 
      dq1=IdealGas2D::VariableDelta<IdealGas2D::VariableType<'0'>>(dq0);

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1.var[0], dq0.var[0], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1.var[1], dq0.var[1], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1.var[2], dq0.var[2], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1.var[3], dq0.var[3], Types::EPS );

      return;
  }

   void Test_IdealGas2D_VariableDelta::test_accessors()
  {
      dq0.var[0]=1.;
      dq0.var[1]=2.;
      dq0.var[2]=3.;
      dq0.var[3]=4.;

   // standard accessor
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0.var[0], dq0[0], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0.var[1], dq0[1], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0.var[2], dq0[2], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0.var[3], dq0[3], Types::EPS );

   // const accessor
      const IdealGas2D::VariableDelta<IdealGas2D::VariableType<'0'>>  dq2(dq0);

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq2[0], dq0[0], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq2[1], dq0[1], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq2[2], dq0[2], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq2[3], dq0[3], Types::EPS );

      return;
  }

   void Test_IdealGas2D_VariableDelta::test_operator_assignment()
  {
      dq0=3.4;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[0], 3.4, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[1], 3.4, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[2], 3.4, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[3], 3.4, Types::EPS );
      return;
  }

   void Test_IdealGas2D_VariableDelta::test_inplace_addition()
  {
      Types::Real  a=2.3;
      dq0=a;

      dq1[0]=-1.;
      dq1[1]= 2.;
      dq1[2]=-3.;
      dq1[3]= 4.;

      dq0+=dq1;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[0], a+dq1[0], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[1], a+dq1[1], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[2], a+dq1[2], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[3], a+dq1[3], Types::EPS );
  }

   void Test_IdealGas2D_VariableDelta::test_inplace_subtraction()
  {
      Types::Real  a=2.3;
      dq0=a;

      dq1[0]=-1.;
      dq1[1]= 2.;
      dq1[2]=-3.;
      dq1[3]= 4.;

      dq0-=dq1;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[0], a-dq1[0], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[1], a-dq1[1], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[2], a-dq1[2], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq0[3], a-dq1[3], Types::EPS );
  }

   void Test_IdealGas2D_VariableDelta::test_inplace_multiplication()
  {
      Types::Real  a=2.3;

      dq1[0]=-1.;
      dq1[1]= 2.;
      dq1[2]=-3.;
      dq1[3]= 4.;

      dq0 =dq1;
      dq1*=a;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1[0], dq0[0]*a, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1[1], dq0[1]*a, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1[2], dq0[2]*a, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1[3], dq0[3]*a, Types::EPS );
  }

   void Test_IdealGas2D_VariableDelta::test_inplace_division()
  {
      Types::Real  a=2.3;

      dq1[0]=-1.;
      dq1[1]= 2.;
      dq1[2]=-3.;
      dq1[3]= 4.;

      dq0 =dq1;
      dq1/=a;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1[0], dq0[0]/a, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1[1], dq0[1]/a, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1[2], dq0[2]/a, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dq1[3], dq0[3]/a, Types::EPS );
  }

