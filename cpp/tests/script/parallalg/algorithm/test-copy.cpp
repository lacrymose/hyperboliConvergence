# include <cppunit/ui/text/TestRunner.h>
# include <cppunit/TestResult.h>

# include <parallalg/algorithm/test-copy.h>

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      runner.addTest( Test_par_copy::suite() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
