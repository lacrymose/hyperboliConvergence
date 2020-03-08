
# include <idealGas2D/test-VariableSet.h>

   void Test_IdealGas2D_VariableSet::setUp(){}

   void Test_IdealGas2D_VariableSet::tearDown(){}

   void Test_IdealGas2D_VariableSet::test_constructors()
  {
   // default constructor
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0.var[0], 0, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0.var[1], 0, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0.var[2], 0, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0.var[3], 0, Types::EPS );

      q0.var[0]=1.;
      q0.var[1]=2.;
      q0.var[2]=3.;
      q0.var[3]=4.;

      q1=IdealGas2D::VariableSet<IdealGas2D::VariableType<'0'>>(q0);

   // copy constructor
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q1.var[0], q0.var[0], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q1.var[1], q0.var[1], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q1.var[2], q0.var[2], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q1.var[3], q0.var[3], Types::EPS );

      return;
  }

   void Test_IdealGas2D_VariableSet::test_accessors()
  {
      q0.var[0]=1.;
      q0.var[1]=2.;
      q0.var[2]=3.;
      q0.var[3]=4.;

   // standard accessor
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0.var[0], q0[0], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0.var[1], q0[1], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0.var[2], q0[2], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0.var[3], q0[3], Types::EPS );

      const IdealGas2D::VariableSet<IdealGas2D::VariableType<'0'>>  q2(q0);

   // const accessor
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q2[0], q0[0], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q2[1], q0[1], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q2[2], q0[2], Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q2[3], q0[3], Types::EPS );
      
      return;
  }

   void Test_IdealGas2D_VariableSet::test_operator_equality()
  {
      q0=3.4;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0[0], 3.4, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0[1], 3.4, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0[2], 3.4, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( q0[3], 3.4, Types::EPS );
      return;
  }

