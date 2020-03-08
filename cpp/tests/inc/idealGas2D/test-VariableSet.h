# ifndef  TEST_IDEALGAS2D_VARIABLESET_H
# define  TEST_IDEALGAS2D_VARIABLESET_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <idealGas2D/idealGas2D.h>

/*
   Tests for public methods of the IdealGas2D_VariableSet class
*/

   class Test_IdealGas2D_VariableSet : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_IdealGas2D_VariableSet );

         CPPUNIT_TEST( test_constructors        );
         CPPUNIT_TEST( test_accessors           );
         CPPUNIT_TEST( test_operator_assignment );

      CPPUNIT_TEST_SUITE_END();

      IdealGas2D::VariableSet<IdealGas2D::VariableType<'0'>>  q0,q1;

   public:
      void setUp();
      void tearDown();

      void test_constructors();
      void test_accessors();

      void test_operator_assignment();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_IdealGas2D_VariableSet );

# endif
