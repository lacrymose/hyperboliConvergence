
# include <array>
# include <iostream>

   template<int N>
   struct Array
  {
      std::array<size_t,N> x;
  };


   int main()
  {
      std::array arr0{1,2,3};

      Array arr1{0,1,2};

      return 0;
  }
