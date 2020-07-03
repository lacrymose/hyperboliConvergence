
# pragma once

# include <cppunit/TestFixture.h>
# include <cppunit/extensions/HelperMacros.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

/*
   Tests copy function of parallalg library
*/

   class Test_par_copy : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_par_copy );

//       CPPUNIT_TEST( test_trivially_copyable );
//       CPPUNIT_TEST( test_not_trivially_copyable );

      CPPUNIT_TEST_SUITE_END();

   public:
      void setUp();
      void tearDown();

//    void test_trivially_copyable();
//    void test_not_trivially_copyable();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_par_copy );

