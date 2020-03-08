# include <cppunit/ui/text/TestRunner.h>
# include <cppunit/TestResult.h>

# include <idealGas2D/test-VariableSet.h>

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      runner.addTest( Test_IdealGas2D_VariableSet::suite() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
