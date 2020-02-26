
# include <iostream>

# include "inc/utils.h"

   int main()
  {
      std::cout << Utils::sign( 1) << std::endl;
      std::cout << Utils::sign(-1) << std::endl;
      std::cout << Utils::sign( 0) << std::endl;

      std::cout << std::endl;

      std::cout << Utils::greaterThan0( 1) << std::endl;
      std::cout << Utils::greaterThan0(-1) << std::endl;
      std::cout << Utils::greaterThan0( 0) << std::endl;

      std::cout << std::endl;

      std::cout << Utils::lessThan0( 1) << std::endl;
      std::cout << Utils::lessThan0(-1) << std::endl;
      std::cout << Utils::lessThan0( 0) << std::endl;

      std::cout << std::endl;

      std::cout << Utils::greaterThan( 4,2) << std::endl;
      std::cout << Utils::greaterThan(-4,2) << std::endl;
      std::cout << Utils::greaterThan( 0,2) << std::endl;

      std::cout << std::endl;

      std::cout << Utils::lessThan( 4,2) << std::endl;
      std::cout << Utils::lessThan(-4,2) << std::endl;
      std::cout << Utils::lessThan( 0,2) << std::endl;

      return 0;
  }
