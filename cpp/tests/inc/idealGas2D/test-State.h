# ifndef  TEST_IDEALGAS2D_SPECIES_H
# define  TEST_IDEALGAS2D_SPECIES_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <idealGas2D/idealGas2D.h>

/*
   Tests for public methods of the IdealGas2D::State class
*/

   class Test_IdealGas2D_State : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_IdealGas2D_State );

         CPPUNIT_TEST( test_accessors );
         CPPUNIT_TEST( test_equality  );

         CPPUNIT_TEST( test_initialisation_conserved );
         CPPUNIT_TEST( test_initialisation_viscous   );

      CPPUNIT_TEST_SUITE_END();

      IdealGas2D::State  s0,s1;

   public:
      void setUp();
      void tearDown();

      void test_accessors();

      void test_equality();

      void test_initialisation_conserved();
      void test_initialisation_viscous();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_IdealGas2D_State );

# endif
