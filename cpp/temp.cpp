
# include <iostream>

   enum struct Numbers : size_t 
  {
      One,
      Two
  };

   int main()
  {
      std::cout << size_t(Numbers::One) << std::endl;

      Numbers three=Numbers(size_t(3));
      Numbers four=Numbers(4);

      std::cout << size_t(three) << std::endl;
      std::cout << size_t(four)  << std::endl;

      size_t four_s = size_t(four);

      std::cout << four_s << std::endl;

      return 0;
  }
